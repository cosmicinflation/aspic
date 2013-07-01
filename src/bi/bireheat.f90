!brane inflation reheating functions in the slow-roll approximations

module bireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf,ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use bisr, only : bi_epsilon_one, bi_epsilon_two, bi_epsilon_three
  use bisr, only : bi_norm_potential
  use bisr, only : bi_x_endinf, bi_efold_primitive
  implicit none

  private

  public bi_x_star, bi_lnrhoreh_max
  public bi_x_rrad, bi_x_rreh


contains

!returns x given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the correspoding bfoldstar
  function bi_x_star(p,mu,w,lnRhoReh,Pstar,bfold)    
    implicit none
    real(kp) :: bi_x_star
    real(kp), intent(in) :: p,mu,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: biData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = bi_x_endinf(p,mu)
    epsOneEnd = bi_epsilon_one(xEnd,p,mu)

    potEnd = bi_norm_potential(xEnd,p,mu)
    primEnd = bi_efold_primitive(xEnd,p,mu)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    biData%real1 = p
    biData%real2 = mu
    biData%real3 = w
    biData%real4 = calF + primEnd

    mini = xend
    maxi = 1._kp/epsilon(1._kp)

    x = zbrent(find_bi_x_star,mini,maxi,tolFind,biData)
    bi_x_star = x

    if (present(bfold)) then
       bfold = -(bi_efold_primitive(x,p,mu) - primEnd)
    endif

    if (x.gt.1._kp) then
       if (display) write(*,*) 'bi_x_star: phi>mu!'
    endif

  end function bi_x_star

  function find_bi_x_star(x,biData)   
    implicit none
    real(kp) :: find_bi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: biData

    real(kp) :: primStar,p,mu,w,CalFplusPrimEnd,potStar,epsOneStar

    p=biData%real1
    mu = biData%real2
    w = biData%real3
    CalFplusPrimEnd = biData%real4

    primStar = bi_efold_primitive(x,p,mu)
    epsOneStar = bi_epsilon_one(x,p,mu)
    potStar = bi_norm_potential(x,p,mu)


    find_bi_x_star = find_reheat(PrimStar,calFplusPrimEnd,w,epsOneStar,potStar)

  end function find_bi_x_star

!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function bi_x_rrad(p,mu,lnRrad,Pstar,bfold)    
    implicit none
    real(kp) :: bi_x_rrad
    real(kp), intent(in) :: p,mu,lnRrad,Pstar
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: biData
    
    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif
   
    xEnd = bi_x_endinf(p,mu)
    epsOneEnd = bi_epsilon_one(xEnd,p,mu)

    potEnd = bi_norm_potential(xEnd,p,mu)
    primEnd = bi_efold_primitive(xEnd,p,mu)
   
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    biData%real1 = p
    biData%real2 = mu
    biData%real3 = calF + primEnd

    mini = xend
    maxi = 1._kp/epsilon(1._kp)

    x = zbrent(find_bi_x_rrad,mini,maxi,tolFind,biData)
    bi_x_rrad = x

    if (present(bfold)) then
       bfold = -(bi_efold_primitive(x,p,mu) - primEnd)
    endif

    if (x.gt.1._kp) then
       if (display) write(*,*) 'bi_x_rrad: phi>mu!'
    endif

  end function bi_x_rrad

  function find_bi_x_rrad(x,biData)   
    implicit none
    real(kp) :: find_bi_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: biData

    real(kp) :: primStar,p,mu,CalFplusPrimEnd,potStar,epsOneStar

    p=biData%real1
    mu = biData%real2
    CalFplusPrimEnd = biData%real3

    primStar = bi_efold_primitive(x,p,mu)
    epsOneStar = bi_epsilon_one(x,p,mu)
    potStar = bi_norm_potential(x,p,mu)

    find_bi_x_rrad = find_reheat_rrad(PrimStar,calFplusPrimEnd,epsOneStar,potStar)

  end function find_bi_x_rrad



!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function bi_x_rreh(p,mu,lnRreh,bfold)    
    implicit none
    real(kp) :: bi_x_rreh
    real(kp), intent(in) :: p,mu,lnRreh
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: biData
    
    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif
   
    xEnd = bi_x_endinf(p,mu)
    epsOneEnd = bi_epsilon_one(xEnd,p,mu)

    potEnd = bi_norm_potential(xEnd,p,mu)
    primEnd = bi_efold_primitive(xEnd,p,mu)
   
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    biData%real1 = p
    biData%real2 = mu
    biData%real3 = calF + primEnd

    mini = xend
    maxi = 1._kp/epsilon(1._kp)

    x = zbrent(find_bi_x_rreh,mini,maxi,tolFind,biData)
    bi_x_rreh = x

    if (present(bfold)) then
       bfold = -(bi_efold_primitive(x,p,mu) - primEnd)
    endif

    if (x.gt.1._kp) then
       if (display) write(*,*) 'bi_x_rreh: phi>mu!'
    endif

  end function bi_x_rreh

  function find_bi_x_rreh(x,biData)   
    implicit none
    real(kp) :: find_bi_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: biData

    real(kp) :: primStar,p,mu,CalFplusPrimEnd,potStar

    p=biData%real1
    mu = biData%real2
    CalFplusPrimEnd = biData%real3

    primStar = bi_efold_primitive(x,p,mu)
    potStar = bi_norm_potential(x,p,mu)

    find_bi_x_rreh = find_reheat_rreh(PrimStar,calFplusPrimEnd,potStar)

  end function find_bi_x_rreh



  function bi_lnrhoreh_max(p,mu,Pstar) 
    implicit none
    real(kp) :: bi_lnrhoreh_max
    real(kp), intent(in) :: p,mu,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp), parameter :: wrad = 1._kp/3._kp
    real(kp), parameter :: junk= 0._kp
    real(kp) :: lnRhoEnd
    
    xEnd = bi_x_endinf(p,mu)       
    potEnd  = bi_norm_potential(xEnd,p,mu)
    epsOneEnd = bi_epsilon_one(xEnd,p,mu)
       
    x = bi_x_star(p,mu,wrad,junk,Pstar)    
    potStar = bi_norm_potential(x,p,mu)
    epsOneStar = bi_epsilon_one(x,p,mu)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'bi_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    bi_lnrhoreh_max = lnRhoEnd

  end function bi_lnrhoreh_max

  
 
end module bireheat

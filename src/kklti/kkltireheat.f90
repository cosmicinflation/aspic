!brane inflation reheating functions in the slow-roll approximations

module kkltireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf,ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use kkltisr, only : kklti_epsilon_one, kklti_epsilon_two, kklti_epsilon_three
  use kkltisr, only : kklti_norm_potential
  use kkltisr, only : kklti_efold_primitive
  implicit none

  private

  public kklti_x_star, kklti_lnrhoreh_max
  public kklti_x_rrad, kklti_x_rreh


contains

!returns x given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the correspoding bfoldstar
  function kklti_x_star(p,mu,xend,w,lnRhoReh,Pstar,bfold)    
    implicit none
    real(kp) :: kklti_x_star
    real(kp), intent(in) :: p,mu,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: kkltiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = kklti_epsilon_one(xEnd,p,mu)

    potEnd = kklti_norm_potential(xEnd,p,mu)
    primEnd = kklti_efold_primitive(xEnd,p,mu)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    kkltiData%real1 = p
    kkltiData%real2 = mu
    kkltiData%real3 = w
    kkltiData%real4 = calF + primEnd

    mini = xend
    maxi = 1._kp/epsilon(1._kp)

    x = zbrent(find_kklti_x_star,mini,maxi,tolFind,kkltiData)
    kklti_x_star = x

    if (present(bfold)) then
       bfold = -(kklti_efold_primitive(x,p,mu) - primEnd)
    endif

    if (x.gt.1._kp) then
       if (display) write(*,*) 'kklti_x_star: phi>mu!'
    endif

  end function kklti_x_star

  function find_kklti_x_star(x,kkltiData)   
    implicit none
    real(kp) :: find_kklti_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: kkltiData

    real(kp) :: primStar,p,mu,w,CalFplusPrimEnd,potStar,epsOneStar

    p=kkltiData%real1
    mu = kkltiData%real2
    w = kkltiData%real3
    CalFplusPrimEnd = kkltiData%real4

    primStar = kklti_efold_primitive(x,p,mu)
    epsOneStar = kklti_epsilon_one(x,p,mu)
    potStar = kklti_norm_potential(x,p,mu)


    find_kklti_x_star = find_reheat(PrimStar,calFplusPrimEnd,w,epsOneStar,potStar)

  end function find_kklti_x_star

!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function kklti_x_rrad(p,mu,xend,lnRrad,Pstar,bfold)    
    implicit none
    real(kp) :: kklti_x_rrad
    real(kp), intent(in) :: p,mu,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: kkltiData
    
    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif
   
    epsOneEnd = kklti_epsilon_one(xEnd,p,mu)

    potEnd = kklti_norm_potential(xEnd,p,mu)
    primEnd = kklti_efold_primitive(xEnd,p,mu)
   
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    kkltiData%real1 = p
    kkltiData%real2 = mu
    kkltiData%real3 = calF + primEnd

    mini = xend
    maxi = 1._kp/epsilon(1._kp)

    x = zbrent(find_kklti_x_rrad,mini,maxi,tolFind,kkltiData)
    kklti_x_rrad = x

    if (present(bfold)) then
       bfold = -(kklti_efold_primitive(x,p,mu) - primEnd)
    endif

    if (x.gt.1._kp) then
       if (display) write(*,*) 'kklti_x_rrad: phi>mu!'
    endif

  end function kklti_x_rrad

  function find_kklti_x_rrad(x,kkltiData)   
    implicit none
    real(kp) :: find_kklti_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: kkltiData

    real(kp) :: primStar,p,mu,CalFplusPrimEnd,potStar,epsOneStar

    p=kkltiData%real1
    mu = kkltiData%real2
    CalFplusPrimEnd = kkltiData%real3

    primStar = kklti_efold_primitive(x,p,mu)
    epsOneStar = kklti_epsilon_one(x,p,mu)
    potStar = kklti_norm_potential(x,p,mu)

    find_kklti_x_rrad = find_reheat_rrad(PrimStar,calFplusPrimEnd,epsOneStar,potStar)

  end function find_kklti_x_rrad



!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function kklti_x_rreh(p,mu,xend,lnRreh,bfold)    
    implicit none
    real(kp) :: kklti_x_rreh
    real(kp), intent(in) :: p,mu,xend,lnRreh
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: kkltiData
    
    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif
   
    epsOneEnd = kklti_epsilon_one(xEnd,p,mu)

    potEnd = kklti_norm_potential(xEnd,p,mu)
    primEnd = kklti_efold_primitive(xEnd,p,mu)
   
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    kkltiData%real1 = p
    kkltiData%real2 = mu
    kkltiData%real3 = calF + primEnd

    mini = xend
    maxi = 1._kp/epsilon(1._kp)

    x = zbrent(find_kklti_x_rreh,mini,maxi,tolFind,kkltiData)
    kklti_x_rreh = x

    if (present(bfold)) then
       bfold = -(kklti_efold_primitive(x,p,mu) - primEnd)
    endif

    if (x.gt.1._kp) then
       if (display) write(*,*) 'kklti_x_rreh: phi>mu!'
    endif

  end function kklti_x_rreh

  function find_kklti_x_rreh(x,kkltiData)   
    implicit none
    real(kp) :: find_kklti_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: kkltiData

    real(kp) :: primStar,p,mu,CalFplusPrimEnd,potStar

    p=kkltiData%real1
    mu = kkltiData%real2
    CalFplusPrimEnd = kkltiData%real3

    primStar = kklti_efold_primitive(x,p,mu)
    potStar = kklti_norm_potential(x,p,mu)

    find_kklti_x_rreh = find_reheat_rreh(PrimStar,calFplusPrimEnd,potStar)

  end function find_kklti_x_rreh



  function kklti_lnrhoreh_max(p,mu,xend,Pstar) 
    implicit none
    real(kp) :: kklti_lnrhoreh_max
    real(kp), intent(in) :: p,mu,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp), parameter :: wrad = 1._kp/3._kp
    real(kp), parameter :: junk= 0._kp
    real(kp) :: lnRhoEnd
    
    potEnd  = kklti_norm_potential(xEnd,p,mu)
    epsOneEnd = kklti_epsilon_one(xEnd,p,mu)
       
    x = kklti_x_star(p,mu,xend,wrad,junk,Pstar)    
    potStar = kklti_norm_potential(x,p,mu)
    epsOneStar = kklti_epsilon_one(x,p,mu)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'kklti_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    kklti_lnrhoreh_max = lnRhoEnd

  end function kklti_lnrhoreh_max

  
 
end module kkltireheat

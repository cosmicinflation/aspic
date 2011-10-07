!running mass reheating functions in the slow-roll

module hyreheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_energy_endinf
  use hysrevol, only : hy_norm_potential, hy_check_slowroll
  use hysrevol, only : hy_epsilon_one, hy_epsilon_two
  use hysrevol, only : hy_x_endinf, hy_nufunc
  implicit none

  private

  public hy_x_reheat, hy_lnrhoend
  public find_hyreheat

contains

!returns x given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the correspoding bfoldstar
  function hy_x_reheat(p,mu,xstop,w,lnRhoReh,Pstar,bfold)    
    implicit none
    real(kp) :: hy_x_reheat
    real(kp), intent(in) :: p,mu,xstop,lnRhoReh,w,Pstar
    real(kp), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: nuEnd,epsEnd,xend,potEnd

    type(transfert) :: hyData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xend = hy_x_endinf(p,mu,xstop)

    epsEnd = hy_epsilon_one(xEnd,p,mu)

    potEnd = hy_norm_potential(xEnd,p)
    nuEnd = hy_nufunc(xEnd,p,mu)
    
!cobe normalised
!    Pstar = quadrupole_to_primscalar(QhysOverT)

    calF = get_calfconst(lnRhoReh,Pstar,w,epsEnd,potEnd)

    hyData%real1 = p
    hyData%real2 = mu
    hyData%real3 = w
    hyData%real4 = calF + nuEnd

    mini = xend
    maxi = 1._kp/epsilon(1._kp)
       
    x = zbrent(find_hyreheat,mini,maxi,tolFind,hyData)

    call hy_check_slowroll(x,xend,p,mu)

    hy_x_reheat = x
    
    if (present(bfold)) then
       bfold = -(hy_nufunc(x,p,mu) - nuEnd)
    endif
    
  end function hy_x_reheat

  function find_hyreheat(x,hyData)   
    implicit none
    real(kp) :: find_hyreheat
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: hyData

    real(kp) :: nuStar,p,mu,w,CalFplusNuEnd,potStar,epsStar

    p=hyData%real1
    mu = hyData%real2
    w = hyData%real3
    CalFplusNuEnd = hyData%real4

    nuStar = hy_nufunc(x,p,mu)
    epsStar = hy_epsilon_one(x,p,mu)
    potStar = hy_norm_potential(x,p)

    find_hyreheat = find_reheat(nuStar,calFplusNuEnd,w,epsStar,potStar)

  end function find_hyreheat



  function hy_lnrhoend(p,mu,xstop,Pstar) 
    implicit none
    real(kp) :: hy_lnrhoend
    real(kp), intent(in) :: p,mu,xstop,Pstar

    real(kp) :: xEnd, potEnd, epsEnd
    real(kp) :: x, potStar, epsStar

    real(kp), parameter :: w = 1._kp/3._kp
    real(kp), parameter :: lnRhoReh = 0._kp
    real(kp) :: lnRhoEnd
    
    xEnd = hy_x_endinf(p,mu,xstop)  
    potEnd  = hy_norm_potential(xEnd,p)
    epsEnd = hy_epsilon_one(xEnd,p,mu)
    
    x = hy_x_reheat(p,mu,xstop,w,lnRhoReh,Pstar)
    potStar = hy_norm_potential(x,p)
    epsStar = hy_epsilon_one(x,p,mu)
    
    if (.not.slowroll_validity(epsStar)) stop 'hy_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_energy_endinf(Pstar,epsStar,epsEnd,potEnd/potStar)

    hy_lnrhoend = lnRhoEnd

  end function hy_lnrhoend

    
end module hyreheat

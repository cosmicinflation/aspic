!KKLT reheating functions in the slow-roll + minimal kinetic term
!approximations

module kkltreheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_energy_endinf
  use kkltsrevol, only : kklt_epsilon_one, kklt_epsilon_two, kklt_norm_potential
  use kkltsrevol, only : kklt_x_endinf, kklt_nufunc
  implicit none

  private

  public kklt_x_star, kklt_lnrhoend
  public find_kklt_x_star

contains

!returns x given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the correspoding bfoldstar
  function kklt_x_star(p,mu,w,lnRhoReh,Pstar,bfold)    
    implicit none
    real(kp) :: kklt_x_star
    real(kp), intent(in) :: p,mu,lnRhoReh,w,Pstar
    real(kp), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: nuEnd,epsEnd,xend,potEnd

    type(transfert) :: kkltData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = kklt_x_endinf(p,mu)
    epsEnd = kklt_epsilon_one(xEnd,p,mu)

    potEnd = kklt_norm_potential(xEnd,p)
    nuEnd = kklt_nufunc(xEnd,p,mu)
   
!cobe normalised
!    Pstar = quadrupole_to_primscalar(QrmsOverT)

    calF = get_calfconst(lnRhoReh,Pstar,w,epsEnd,potEnd)

    kkltData%real1 = p
    kkltData%real2 = mu
    kkltData%real3 = w
    kkltData%real4 = calF + nuEnd

    mini = xend
    maxi = 1._kp/epsilon(1._kp)

    x = zbrent(find_kklt_x_star,mini,maxi,tolFind,kkltData)
    kklt_x_star = x

    if (present(bfold)) then
       bfold = -(kklt_nufunc(x,p,mu) - nuEnd)
    endif
    
  end function kklt_x_star

  function find_kklt_x_star(x,kkltData)   
    implicit none
    real(kp) :: find_kklt_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: kkltData

    real(kp) :: nuStar,p,mu,w,CalFplusNuEnd,potStar,epsStar

    p=kkltData%real1
    mu = kkltData%real2
    w = kkltData%real3
    CalFplusNuEnd = kkltData%real4

    nuStar = kklt_nufunc(x,p,mu)
    epsStar = kklt_epsilon_one(x,p,mu)
    potStar = kklt_norm_potential(x,p)

    find_kklt_x_star = find_reheat(nuStar,calFplusNuEnd,w,epsStar,potStar)

  end function find_kklt_x_star



  function kklt_lnrhoend(p,mu,Pstar) 
    implicit none
    real(kp) :: kklt_lnrhoend
    real(kp), intent(in) :: p,mu,Pstar

    real(kp) :: xEnd, potEnd, epsEnd
    real(kp) :: x, potStar, epsStar

    real(kp), parameter :: w = 1._kp/3._kp
    real(kp), parameter :: lnRhoReh = 0._kp
    real(kp) :: lnRhoEnd
    
    xEnd = kklt_x_endinf(p,mu)      
    potEnd  = kklt_norm_potential(xEnd,p)
    epsEnd = kklt_epsilon_one(xEnd,p,mu)
       
    x = kklt_x_star(p,mu,w,lnRhoReh,Pstar)
    potStar = kklt_norm_potential(x,p)
    epsStar = kklt_epsilon_one(x,p,mu)
    
    if (.not.slowroll_validity(epsStar)) stop 'kklt_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_energy_endinf(Pstar,epsStar,epsEnd,potEnd/potStar)

    kklt_lnrhoend = lnRhoEnd

  end function kklt_lnrhoend

    
end module kkltreheat

!kklt as a small field model: reheating functions in
!the slow-roll approximations

module kksfreheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent  
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_energy_endinf
  use kksfsrevol, only : kksf_epsilon_one, kksf_epsilon_two, kksf_norm_potential
  use kksfsrevol, only : kksf_x_endinf, kksf_nufunc
  use sfreheat, only : find_sfi_x_star
  implicit none

  private

  public kksf_x_star, kksf_lnrhoend

contains

!returns x given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the correspoding bfoldstar
  function kksf_x_star(p,mu,w,lnRhoReh,Pstar,bfold)    
    implicit none
    real(kp) :: kksf_x_star
    real(kp), intent(in) :: p,mu,lnRhoReh,w,Pstar
    real(kp), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: nuEnd,epsEnd,xend,potEnd

    type(transfert) :: kksfData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = kksf_x_endinf(p,mu)
    epsEnd = kksf_epsilon_one(xEnd,p,mu)

    potEnd = kksf_norm_potential(xEnd,p)
    nuEnd = kksf_nufunc(xEnd,p,mu)
   
!cobe normalised
!    Pstar = quadrupole_to_primscalar(QrmsOverT)

    calF = get_calfconst(lnRhoReh,Pstar,w,epsEnd,potEnd)

    kksfData%real1 = -p
    kksfData%real2 = mu
    kksfData%real3 = w
    kksfData%real4 = calF + nuEnd

    mini = xEnd
    maxi = 1._kp/epsilon(1._kp)

    x = zbrent(find_sfi_x_star,mini,maxi,tolFind,kksfData)
    kksf_x_star = x

    if (present(bfold)) then
       bfold = -(kksf_nufunc(x,p,mu) - nuEnd)
    endif

    if (x.lt.1._kp) then
       if (display) write(*,*) 'kksf_x_star: phi<mu!'
    endif
  end function kksf_x_star



  function kksf_lnrhoend(p,mu,Pstar) 
    implicit none
    real(kp) :: kksf_lnrhoend
    real(kp), intent(in) :: p,mu,Pstar

    real(kp) :: xEnd, potEnd, epsEnd
    real(kp) :: x, potStar, epsStar

    real(kp), parameter :: w = 1._kp/3._kp
    real(kp), parameter :: lnRhoReh = 0._kp
    real(kp) :: lnRhoEnd
    
    xEnd = kksf_x_endinf(p,mu)       
    potEnd  = kksf_norm_potential(xEnd,p)
    epsEnd = kksf_epsilon_one(xEnd,p,mu)
       
    x = kksf_x_star(p,mu,w,lnRhoReh,Pstar)    
    potStar = kksf_norm_potential(x,p)
    epsStar = kksf_epsilon_one(x,p,mu)
    
    if (.not.slowroll_validity(epsStar)) stop 'kksf_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_energy_endinf(Pstar,epsStar,epsEnd,potEnd/potStar)

    kksf_lnrhoend = lnRhoEnd

  end function kksf_lnrhoend

    

end module kksfreheat

!running mass reheating functions in the slow-roll

module rmreheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_energy_endinf
  use rmsrevol, only : rm_norm_potential
  use rmsrevol, only : rm_epsilon_one, rm_epsilon_two
  use rmsrevol, only : rm_x_endinf, rm_nufunc
  implicit none

  private

  public rm_x_reheat, rm_lnrhoend
  public find_rmreheat

contains

!returns x given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the correspoding bfoldstar
  function rm_x_reheat(p,mu,nu,xstop,w,lnRhoReh,Pstar,bfold)    
    implicit none
    real(kp) :: rm_x_reheat
    real(kp), intent(in) :: p,mu,nu,xstop,lnRhoReh,w,Pstar
    real(kp), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: nuEnd,epsEnd,xend,potEnd

    type(transfert) :: rmData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xend = rm_x_endinf(p,mu,nu,xstop)

    epsEnd = rm_epsilon_one(xEnd,p,mu,nu)

    potEnd = rm_norm_potential(xEnd,p,mu,nu)
    nuEnd = rm_nufunc(xEnd,p,mu,nu)
   
!cobe normalised
!    Pstar = quadrupole_to_primscalar(QrmsOverT)

    calF = get_calfconst(lnRhoReh,Pstar,w,epsEnd,potEnd)

    rmData%real1 = p
    rmData%real2 = mu
    rmData%real3 = nu
    rmData%real4 = w
    rmData%real5 = calF + nuEnd

    if (xend.lt.1._kp) then
       mini = epsilon(1._kp)
       maxi = 1._kp - epsilon(1._kp)
    elseif (xend.gt.1._kp) then
       mini = 1._kp + epsilon(1._kp)
       if (nu.gt.0._kp) then
          maxi = xend
       else
          maxi = 1._kp/epsilon(1._kp)
       endif
    else
       stop 'rm_x_reheat: internal error'
    endif
    
    x = zbrent(find_rmreheat,mini,maxi,tolFind,rmData)

    rm_x_reheat = x

    if (present(bfold)) then
       bfold = -(rm_nufunc(x,p,mu,nu) - nuEnd)
    endif
    
  end function rm_x_reheat

  function find_rmreheat(x,rmData)   
    implicit none
    real(kp) :: find_rmreheat
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rmData

    real(kp) :: nuStar,p,mu,nu,w,CalFplusNuEnd,potStar,epsStar

    p=rmData%real1
    mu = rmData%real2
    nu = rmData%real3
    w = rmData%real4
    CalFplusNuEnd = rmData%real5

    nuStar = rm_nufunc(x,p,mu,nu)
    epsStar = rm_epsilon_one(x,p,mu,nu)
    potStar = rm_norm_potential(x,p,mu,nu)

    find_rmreheat = find_reheat(nuStar,calFplusNuEnd,w,epsStar,potStar)

  end function find_rmreheat



  function rm_lnrhoend(p,mu,nu,xstop,Pstar) 
    implicit none
    real(kp) :: rm_lnrhoend
    real(kp), intent(in) :: p,mu,nu,xstop,Pstar

    real(kp) :: xEnd, potEnd, epsEnd
    real(kp) :: x, potStar, epsStar

    real(kp), parameter :: w = 1._kp/3._kp
    real(kp), parameter :: lnRhoReh = 0._kp
    real(kp) :: lnRhoEnd
    
    xEnd = rm_x_endinf(p,mu,nu,xstop)  
    potEnd  = rm_norm_potential(xEnd,p,mu,nu)
    epsEnd = rm_epsilon_one(xEnd,p,mu,nu)
       
    x = rm_x_reheat(p,mu,nu,xstop,w,lnRhoReh,Pstar)
    potStar = rm_norm_potential(x,p,mu,nu)
    epsStar = rm_epsilon_one(x,p,mu,nu)
    
    if (.not.slowroll_validity(epsStar)) stop 'rm_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_energy_endinf(Pstar,epsStar,epsEnd,potEnd/potStar)

    rm_lnrhoend = lnRhoEnd

  end function rm_lnrhoend

    
end module rmreheat

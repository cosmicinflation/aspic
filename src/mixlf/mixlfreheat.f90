!Mixed large field reheating functions in the slow-roll

module mixlfreheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_energy_endinf
  use mixlfsrevol, only : mixlf_norm_potential
  use mixlfsrevol, only : mixlf_epsilon_one, mixlf_epsilon_two
  use mixlfsrevol, only : mixlf_x_endinf, mixlf_nufunc
  implicit none

  private

  public mixlf_x_reheat, mixlf_lnrhoend
  public find_mixlfreheat

contains

!returns x given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the correspoding bfoldstar
  function mixlf_x_reheat(p,q,alpha,w,lnRhoReh,Pstar,bfold)    
    implicit none
    real(kp) :: mixlf_x_reheat
    real(kp), intent(in) :: p,q,alpha,lnRhoReh,w,Pstar
    real(kp), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: nuEnd,epsEnd,xend,potEnd

    type(transfert) :: mixlfData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = mixlf_x_endinf(p,q,alpha)
    epsEnd = mixlf_epsilon_one(xEnd,p,q,alpha)

    potEnd = mixlf_norm_potential(xEnd,p,q,alpha)
    nuEnd = mixlf_nufunc(xEnd,p,q,alpha)
   
!cobe normalised
!    Pstar = quadrupole_to_primscalar(QrmsOverT)

    calF = get_calfconst(lnRhoReh,Pstar,w,epsEnd,potEnd)

    mixlfData%real1 = p
    mixlfData%real2 = q
    mixlfData%real3 = alpha
    mixlfData%real4 = w
    mixlfData%real5 = calF + nuEnd

    mini = xend
    maxi = 1._kp/epsilon(1._kp)

    x = zbrent(find_mixlfreheat,mini,maxi,tolFind,mixlfData)
    mixlf_x_reheat = x

    if (present(bfold)) then
       bfold = -(mixlf_nufunc(x,p,q,alpha) - nuEnd)
    endif
    
  end function mixlf_x_reheat

  function find_mixlfreheat(x,mixlfData)   
    implicit none
    real(kp) :: find_mixlfreheat
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: mixlfData

    real(kp) :: nuStar,p,q,alpha,w,CalFplusNuEnd,potStar,epsStar

    p=mixlfData%real1
    q = mixlfData%real2
    alpha = mixlfData%real3
    w = mixlfData%real4
    CalFplusNuEnd = mixlfData%real5

    nuStar = mixlf_nufunc(x,p,q,alpha)
    epsStar = mixlf_epsilon_one(x,p,q,alpha)
    potStar = mixlf_norm_potential(x,p,q,alpha)

    find_mixlfreheat = find_reheat(nuStar,calFplusNuEnd,w,epsStar,potStar)

  end function find_mixlfreheat



  function mixlf_lnrhoend(p,q,alpha,Pstar) 
    implicit none
    real(kp) :: mixlf_lnrhoend
    real(kp), intent(in) :: p,q,alpha,Pstar

    real(kp) :: xEnd, potEnd, epsEnd
    real(kp) :: x, potStar, epsStar

    real(kp), parameter :: w = 1._kp/3._kp
    real(kp), parameter :: lnRhoReh = 0._kp
    real(kp) :: lnRhoEnd
    
    xEnd = mixlf_x_endinf(p,q,alpha)      
    potEnd  = mixlf_norm_potential(xEnd,p,q,alpha)
    epsEnd = mixlf_epsilon_one(xEnd,p,q,alpha)
       
    x = mixlf_x_reheat(p,q,alpha,w,lnRhoReh,Pstar)
    potStar = mixlf_norm_potential(x,p,q,alpha)
    epsStar = mixlf_epsilon_one(x,p,q,alpha)
    
    if (.not.slowroll_validity(epsStar)) stop 'mixlf_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_energy_endinf(Pstar,epsStar,epsEnd,potEnd/potStar)

    mixlf_lnrhoend = lnRhoEnd

  end function mixlf_lnrhoend

    
end module mixlfreheat

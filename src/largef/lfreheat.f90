!large field model reheating functions in the slow-roll approximations

module lfreheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_energy_endinf
  use lfsrevol, only : lf_epsilon_one, lf_epsilon_two, lf_norm_potential
  use lfsrevol, only : lf_x_endinf, lf_nufunc
  implicit none

  private

  public lf_x_reheat, lf_lnrhoend, lf_lnrhoreh, lf_x_obs

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the correspoding bfoldstar
  function lf_x_reheat(p,w,lnRhoReh,Pstar,bfold)    
    implicit none
    real(kp) :: lf_x_reheat
    real(kp), intent(in) :: p,lnRhoReh,w,Pstar
    real(kp), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: nuEnd,epsEnd,xend,potEnd

    type(transfert) :: lfData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = lf_x_endinf(p)
    epsEnd = lf_epsilon_one(xEnd,p)

    potEnd = lf_norm_potential(xEnd,p)
    nuEnd = lf_nufunc(xEnd,p)
   

    calF = get_calfconst(lnRhoReh,Pstar,w,epsEnd,potEnd)

    lfData%real1 = p    
    lfData%real2 = w
    lfData%real3 = calF + nuEnd

    mini = xEnd
    maxi = 1._kp/epsilon(1._kp)

    x = zbrent(find_lfreheat,mini,maxi,tolFind,lfData)
    lf_x_reheat = x

    if (present(bfold)) then
       bfold = - (lf_nufunc(x,p) - nuEnd)
    endif

  end function lf_x_reheat

  function find_lfreheat(x,lfData)   
    implicit none
    real(kp) :: find_lfreheat
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: lfData

    real(kp) :: nuStar,p,w,CalFplusNuEnd,potStar,epsStar

    p=lfData%real1
    w = lfData%real2
    CalFplusNuEnd = lfData%real3

    nuStar = lf_nufunc(x,p)
    epsStar = lf_epsilon_one(x,p)
    potStar = lf_norm_potential(x,p)

    find_lfreheat = find_reheat(nuStar,calFplusNuEnd,w,epsStar,potStar)
  end function find_lfreheat



  function lf_lnrhoend(p,Pstar) 
    implicit none
    real(kp) :: lf_lnrhoend
    real(kp), intent(in) :: p,Pstar

    real(kp) :: xEnd, potEnd, epsEnd
    real(kp) :: x, potStar, epsStar

    real(kp), parameter :: w = 1._kp/3._kp
    real(kp), parameter :: lnRhoReh = 0._kp
    real(kp) :: lnRhoEnd
    
    xEnd = lf_x_endinf(p)
    potEnd  = lf_norm_potential(xEnd,p)
    epsEnd = lf_epsilon_one(xEnd,p)
       
    x = lf_x_reheat(p,w,lnRhoReh,Pstar)    
    potStar = lf_norm_potential(x,p)
    epsStar = lf_epsilon_one(x,p)
    
    if (.not.slowroll_validity(epsStar)) stop 'lf_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_energy_endinf(Pstar,epsStar,epsEnd,potEnd/potStar)

    lf_lnrhoend = lnRhoEnd

  end function lf_lnrhoend

  
   
!return the unique p,x giving eps12 (and bfold if input)
  function lf_x_obs(eps1,eps2,bfold)    
    implicit none
    real(kp), dimension(2) :: lf_x_obs
    real(kp), intent(in) :: eps1,eps2
    real(kp), optional :: bfold
    
    real(kp) :: x, xEnd, p

    if (eps2.le.0._kp) then
       stop 'lf_x_obs: eps2<=0'
    endif

    p = 4._kp*eps1/eps2
    x = sqrt(8._kp*eps1/eps2/eps2)

    lf_x_obs(1) = p
    lf_x_obs(2) = x
    
    xEnd = lf_x_endinf(p)    
    
    if (present(bfold)) then        
       bfold = -(lf_nufunc(x,p)-lf_nufunc(xEnd,p))
    endif

  end function lf_x_obs

  
!returns lnrhoreh from eps12, wreh and Pstar (and bfoldstar if
!present)
  function lf_lnrhoreh(wreh,eps1,eps2,Pstar,bfold)     
    implicit none
    real(kp) :: lf_lnrhoreh
    real(kp), intent(inout) :: wreh
    real(kp), intent(in) :: eps1,eps2,Pstar
    real(kp), optional :: bfold

    real(kp) :: w, p
    real(kp) :: xEnd, potEnd, epsEnd
    real(kp) :: x, potStar, epsStar, lnH2OverEps
    real(kp) :: oneMinusThreeWlnTreh, deltaNstar
    real(kp), dimension(2) :: lfstar

    logical, parameter :: printTest = .false.
    logical, parameter :: enforceWofp = .true.

    lfstar = lf_x_obs(eps1,eps2,bfold)    
    p = lfstar(1)
    x = lfstar(2)

    if (enforceWofp) then
       w = (p-2._kp)/(p+2._kp)
       wreh = w
    else
       w = wreh
    endif

    xEnd = lf_x_endinf(p)       
    potEnd  = lf_norm_potential(xEnd,p)
    epsEnd = lf_epsilon_one(xEnd,p)
    
    potStar = lf_norm_potential(x,p)
    epsStar = lf_epsilon_one(x,p)

    if (.not.slowroll_validity(epsStar)) stop 'lf_lnrhoreh: cannot trust slow-roll!'
    
    lnH2OverEps = log(Pstar*8*pi*pi)
           
    deltaNstar = lf_nufunc(x,p) - lf_nufunc(xEnd,p)

    if (printTest) then
       write(*,*)'eps1in= eps1comp= ',eps1, epsStar
       write(*,*)'eps2in= eps2comp= ',eps2,lf_epsilon_two(x,p)
    endif
    
       
    oneMinusThreeWlnTreh &
         = (3._kp + 3._kp*w)*(Nzero + deltaNstar) - 0.5_kp*(1._kp+3._kp*w)*lnH2OverEps &
         + 0.5_kp*(1._kp - 3._kp*w)*log(epsStar) &
         + log(3._kp*(3._kp - epsStar)/(epsStar*(3._kp - epsEnd))) &
         + log(potEnd/potStar)

    lf_lnrhoreh = oneMinusThreeWlnTreh*4._kp/(1._kp - 3._kp*w)
    
  end function lf_lnrhoreh



end module lfreheat

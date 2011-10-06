!small field model reheating functions in the slow-roll approximations

module sfreheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_energy_endinf
  use sfsrevol, only : sf_epsilon_one, sf_epsilon_two, sf_norm_potential
  use sfsrevol, only : sf_x_endinf, sf_nufunc
  implicit none

  private

  public sf_x_reheat, sf_lnrhoend, sf_lnrhoreh, sf_x_obs
  public find_sfreheat

contains

!returns x given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the correspoding bfoldstar
  function sf_x_reheat(p,mu,w,lnRhoReh,Pstar,bfold)    
    implicit none
    real(kp) :: sf_x_reheat
    real(kp), intent(in) :: p,mu,lnRhoReh,w,Pstar
    real(kp), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: nuEnd,epsEnd,xend,potEnd

    type(transfert) :: sfData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = sf_x_endinf(p,mu)
    epsEnd = sf_epsilon_one(xEnd,p,mu)

    potEnd = sf_norm_potential(xEnd,p)
    nuEnd = sf_nufunc(xEnd,p,mu)
   
!cobe normalised
!    Pstar = quadrupole_to_primscalar(QrmsOverT)

    calF = get_calfconst(lnRhoReh,Pstar,w,epsEnd,potEnd)

    sfData%real1 = p
    sfData%real2 = mu
    sfData%real3 = w
    sfData%real4 = calF + nuEnd

    mini = 0._kp
    maxi = xEnd + epsilon(1._kp)

    x = zbrent(find_sfreheat,mini,maxi,tolFind,sfData)
    sf_x_reheat = x

    if (present(bfold)) then
       bfold = -(sf_nufunc(x,p,mu) - nuEnd)
    endif

    if (x.gt.1._kp) then
       if (display) write(*,*) 'sf_x_reheat: phi>mu!'
    endif

  end function sf_x_reheat

  function find_sfreheat(x,sfData)   
    implicit none
    real(kp) :: find_sfreheat
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: sfData

    real(kp) :: nuStar,p,mu,w,CalFplusNuEnd,potStar,epsStar

    p=sfData%real1
    mu = sfData%real2
    w = sfData%real3
    CalFplusNuEnd = sfData%real4

    nuStar = sf_nufunc(x,p,mu)
    epsStar = sf_epsilon_one(x,p,mu)
    potStar = sf_norm_potential(x,p)


    find_sfreheat = find_reheat(nuStar,calFplusNuEnd,w,epsStar,potStar)

!    sr_find_reheat_sf = nuStar - CalFPlusNuCalEnd &
!         + 1._kp/(3._kp+3._kp*w) &
!         * log( (9._kp-3._kp*epsStar)/( 9._kp*(2._kp*epsStar)**(0.5_kp+1.5_kp*w)*potStar) )

  end function find_sfreheat



  function sf_lnrhoend(p,mu,Pstar) 
    implicit none
    real(kp) :: sf_lnrhoend
    real(kp), intent(in) :: p,mu,Pstar

    real(kp) :: xEnd, potEnd, epsEnd
    real(kp) :: x, potStar, epsStar

    real(kp), parameter :: w = 1._kp/3._kp
    real(kp), parameter :: lnRhoReh = 0._kp
    real(kp) :: lnRhoEnd
    
    xEnd = sf_x_endinf(p,mu)       
    potEnd  = sf_norm_potential(xEnd,p)
    epsEnd = sf_epsilon_one(xEnd,p,mu)
       
    x = sf_x_reheat(p,mu,w,lnRhoReh,Pstar)    
    potStar = sf_norm_potential(x,p)
    epsStar = sf_epsilon_one(x,p,mu)
    
    if (.not.slowroll_validity(epsStar)) stop 'sf_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_energy_endinf(Pstar,epsStar,epsEnd,potEnd/potStar)

    sf_lnrhoend = lnRhoEnd

  end function sf_lnrhoend

  
   
!return the unique p,mu,x giving eps123 (and bfold if input)
  function sf_x_obs(eps1,eps2,eps3,bfold)    
    implicit none
    real(kp), dimension(3) :: sf_x_obs
    real(kp), intent(in) :: eps1,eps2,eps3
    real(kp), optional :: bfold

    real(kp), parameter :: tolFind = tolkp   
    type(transfert) :: sfData
    real(kp) :: p,mu
    real(kp) :: x, xToP
    real(kp) :: xEnd

    if ((4._kp*eps1 + eps2 - 2._kp*eps3).eq.0._kp) then      
       stop 'sf_x_obs: p -> infty!'
    endif
    
    p = (2._kp/eps2)*( 8._kp*eps1**2 - 2._kp*eps1*eps2 + eps2**2 - eps2*eps3) &
         /(4._kp*eps1 + eps2 - 2._kp*eps3)

    if (p.lt.0._kp) then
       stop 'sr_x_obs_sf: p < 0!'
    endif


    xToP = 2._kp*eps1*(eps2 - 4._kp*eps1)/(eps2*(eps2 - eps3))

   
    if ((xToP.le.0._kp).or.(xToP.gt.1._kp)) then
       write(*,*)'(phi/mu)^p= ',xToP
       stop 'sf_x_obs: (phi/mu)^p not consistent with small fields'
    endif
    
    x = xToP**(1._kp/p)

    mu = p/sqrt(2._kp*eps1) * x**(p-1)/(1._kp - xToP)
    
    xEnd = sf_x_endinf(p,mu)

    if (present(bfold)) then        
       bfold = -(sf_nufunc(x,p,mu)-sf_nufunc(xEnd,p,mu))
    endif

    sf_x_obs(1) = x
    sf_x_obs(2) = p
    sf_x_obs(3) = mu

  end function sf_x_obs

  !check if eps123 uniquely define x,mu,p
  function sf_checkobs(eps1,eps2,eps3)
    real(kp), intent(in) :: eps1,eps2,eps3
    logical :: sf_checkobs

    real(kp) :: xtop, p

    sf_checkobs = .false.
    
    if (eps2.lt.0._kp) return
    
    xtop = 2._kp*eps1*(eps2 - 4._kp*eps1)/(eps2*(eps2 - eps3))

    if ((xtop.gt.1._kp).or.(xtop.lt.0._kp)) return

    p = (2._kp/eps2)*( 8._kp*eps1**2 - 2._kp*eps1*eps2 + eps2**2 - eps2*eps3) &
         /(4._kp*eps1 + eps2 - 2._kp*eps3)

    if (p.lt.1._kp) return

    sf_checkobs = .true.

  end function sf_checkobs


!returns lnrhoreh from eps123, wreh and Pstar (and bfoldstar if
!present)
  function sf_lnrhoreh(wreh,eps1,eps2,eps3,Pstar,bfold)     
    implicit none
    real(kp) :: sf_lnrhoreh
    real(kp), intent(in) :: wreh,eps1,eps2,eps3,Pstar
    real(kp), optional :: bfold

    real(kp) :: w, p, mu
    real(kp) :: xEnd, potEnd, epsEnd
    real(kp) :: x, potStar, epsStar, lnH2OverEps
    real(kp) :: oneMinusThreeWlnTreh, deltaNstar
    real(kp), dimension(3) :: sfstar

    logical, parameter :: printTest = .false.

    if (wreh.eq.1._kp/3._kp) then
       write(*,*)'sf_lnrhoreh: w = 1/3!'
       stop
    else
       w = wreh
    endif
    
    sfstar = sf_x_obs(eps1,eps2,eps3,bfold)

    x = sfstar(1)
    p = sfstar(2)
    mu = sfstar(3)

    xEnd = sf_x_endinf(p,mu)       
    potEnd  = sf_norm_potential(xEnd,p)
    epsEnd = sf_epsilon_one(xEnd,p,mu)
    
    potStar = sf_norm_potential(x,p)
    epsStar = sf_epsilon_one(x,p,mu)

    if (.not.slowroll_validity(epsStar)) stop 'sf_lnrhoreh: slow-roll violated!'

    lnH2OverEps = log(Pstar*8*pi*pi)
           
    deltaNstar = sf_nufunc(x,p,mu) - sf_nufunc(xEnd,p,mu)

    if (printTest) then
       write(*,*)'eps1in= eps1comp= ',eps1, epsStar
       write(*,*)'eps2in= eps2comp= ',eps2,sf_epsilon_two(x,p,mu)
    endif
    
    oneMinusThreeWlnTreh &
         = (3._kp + 3._kp*w)*(Nzero + deltaNstar) - 0.5_kp*(1._kp+3._kp*w)*lnH2OverEps &
         + 0.5_kp*(1._kp - 3._kp*w)*log(epsStar) &
         + log(3._kp*(3._kp - epsStar)/(epsStar*(3._kp - epsEnd))) &
         + log(potEnd/potStar)

    sf_lnrhoreh = oneMinusThreeWlnTreh*4._kp/(1._kp - 3._kp*w)

  end function sf_lnrhoreh



end module sfreheat

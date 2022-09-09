!hyperbolic model reheating functions in the slow-roll approximations

module hbireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf,ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use hbisr, only : hbi_epsilon_one, hbi_epsilon_two, hbi_epsilon_three
  use hbisr, only : hbi_norm_potential, hbi_numacc_x_potbig
  use hbisr, only : hbi_x_endinf, hbi_efold_primitive
  implicit none

  private

  public hbi_x_star, hbi_lnrhoreh_max, hbi_lnrhoreh_fromepsilon, hbi_xnmu_fromepsilon
  public hbi_x_rrad, hbi_x_rreh


contains


!returns x given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the correspoding bfoldstar
  function hbi_x_star(n,mu,xend,w,lnRhoReh,Pstar,bfold)
    implicit none
    real(kp) :: hbi_x_star
    real(kp), intent(in) :: n,mu,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: hbiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = hbi_epsilon_one(xEnd,n,mu)

    potEnd = hbi_norm_potential(xEnd,n,mu)
    primEnd = hbi_efold_primitive(xEnd,n,mu)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    hbiData%real1 = n
    hbiData%real2 = mu
    hbiData%real3 = w
    hbiData%real4 = calF + primEnd

    mini = xEnd + epsilon(1._kp)
    maxi = hbi_numacc_x_potbig(n)
    
    x = zbrent(find_hbi_x_star,mini,maxi,tolFind,hbiData)
    hbi_x_star = x

    if (present(bfold)) then
       bfold = -(hbi_efold_primitive(x,n,mu) - primEnd)
    endif

    if (x.gt.1._kp) then
       if (display) write(*,*) 'hbi_x_star: phi>mu!'
    endif

  end function hbi_x_star

  function find_hbi_x_star(x,hbiData)   
    implicit none
    real(kp) :: find_hbi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: hbiData

    real(kp) :: primStar,n,mu,w,CalFplusPrimEnd,potStar,epsOneStar

    n=hbiData%real1
    mu = hbiData%real2
    w = hbiData%real3
    CalFplusPrimEnd = hbiData%real4

    primStar = hbi_efold_primitive(x,n,mu)
    epsOneStar = hbi_epsilon_one(x,n,mu)
    potStar = hbi_norm_potential(x,n,mu)
    
    find_hbi_x_star = find_reheat(PrimStar,calFplusPrimEnd,w,epsOneStar,potStar)

  end function find_hbi_x_star

!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function hbi_x_rrad(n,mu,xend,lnRrad,Pstar,bfold)
    implicit none
    real(kp) :: hbi_x_rrad
    real(kp), intent(in) :: n,mu,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: hbiData
    
    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif
   
    epsOneEnd = hbi_epsilon_one(xEnd,n,mu)

    potEnd = hbi_norm_potential(xEnd,n,mu)
    primEnd = hbi_efold_primitive(xEnd,n,mu)
   
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    hbiData%real1 = n
    hbiData%real2 = mu
    hbiData%real3 = calF + primEnd

    mini = xEnd + epsilon(1._kp)
    maxi = hbi_numacc_x_potbig(n)

    x = zbrent(find_hbi_x_rrad,mini,maxi,tolFind,hbiData)
    hbi_x_rrad = x

    if (present(bfold)) then
       bfold = -(hbi_efold_primitive(x,n,mu) - primEnd)
    endif

    if (x.gt.1._kp) then
       if (display) write(*,*) 'hbi_x_rrad: phi>mu!'
    endif

  end function hbi_x_rrad

  function find_hbi_x_rrad(x,hbiData)   
    implicit none
    real(kp) :: find_hbi_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: hbiData

    real(kp) :: primStar,n,mu,CalFplusPrimEnd,potStar,epsOneStar

    n=hbiData%real1
    mu = hbiData%real2
    CalFplusPrimEnd = hbiData%real3

    primStar = hbi_efold_primitive(x,n,mu)
    epsOneStar = hbi_epsilon_one(x,n,mu)
    potStar = hbi_norm_potential(x,n,mu)

    find_hbi_x_rrad = find_reheat_rrad(PrimStar,calFplusPrimEnd,epsOneStar,potStar)

  end function find_hbi_x_rrad



!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function hbi_x_rreh(n,mu,xend,lnRreh,bfold)    
    implicit none
    real(kp) :: hbi_x_rreh
    real(kp), intent(in) :: n,mu,xend,lnRreh
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: hbiData
    
    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif
   
    epsOneEnd = hbi_epsilon_one(xEnd,n,mu)

    potEnd = hbi_norm_potential(xEnd,n,mu)
    primEnd = hbi_efold_primitive(xEnd,n,mu)
   
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    hbiData%real1 = n
    hbiData%real2 = mu
    hbiData%real3 = calF + primEnd

    mini = xEnd
    maxi = hbi_numacc_x_potbig(n)

    x = zbrent(find_hbi_x_rreh,mini,maxi,tolFind,hbiData)
    hbi_x_rreh = x

    if (present(bfold)) then
       bfold = -(hbi_efold_primitive(x,n,mu) - primEnd)
    endif

    if (x.gt.1._kp) then
       if (display) write(*,*) 'hbi_x_rreh: phi>mu!'
    endif

  end function hbi_x_rreh

  function find_hbi_x_rreh(x,hbiData)   
    implicit none
    real(kp) :: find_hbi_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: hbiData

    real(kp) :: primStar,n,mu,CalFplusPrimEnd,potStar

    n=hbiData%real1
    mu = hbiData%real2
    CalFplusPrimEnd = hbiData%real3

    primStar = hbi_efold_primitive(x,n,mu)
    potStar = hbi_norm_potential(x,n,mu)

    find_hbi_x_rreh = find_reheat_rreh(PrimStar,calFplusPrimEnd,potStar)

  end function find_hbi_x_rreh



  function hbi_lnrhoreh_max(n,mu,xend,Pstar) 
    implicit none
    real(kp) :: hbi_lnrhoreh_max
    real(kp), intent(in) :: n,mu,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp), parameter :: wrad = 1._kp/3._kp
    real(kp), parameter :: junk= 0._kp
    real(kp) :: lnRhoEnd
    
    potEnd  = hbi_norm_potential(xEnd,n,mu)
    epsOneEnd = hbi_epsilon_one(xEnd,n,mu)
       
    x = hbi_x_star(n,mu,xend,wrad,junk,Pstar)    
    potStar = hbi_norm_potential(x,n,mu)
    epsOneStar = hbi_epsilon_one(x,n,mu)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'hbi_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    hbi_lnrhoreh_max = lnRhoEnd

  end function hbi_lnrhoreh_max

  
   
!return the unique x,n,mu giving eps123 (and bfold if input)
  function hbi_xnmu_fromepsilon(eps1,eps2,eps3,bfold)
    implicit none
    real(kp), dimension(3) :: hbi_xnmu_fromepsilon
    real(kp), intent(in) :: eps1,eps2,eps3
    real(kp), intent(out), optional :: bfold

    real(kp) :: x, n, mu, xENd

    x = acosh(sqrt(eps3/eps2))

    n = 4._kp*eps1/eps3

    mu = n/sqrt(2._kp*eps1)/tanh(x)

    if (present(bfold)) then        
       xEnd = hbi_x_endinf(n,mu)
       bfold = -(hbi_efold_primitive(x,n,mu)-hbi_efold_primitive(xEnd,n,mu))
    endif

    hbi_xnmu_fromepsilon(1) = x
    hbi_xnmu_fromepsilon(2) = n
    hbi_xnmu_fromepsilon(3) = mu

  end function hbi_xnmu_fromepsilon



!returns lnrhoreh from eps123, wreh and Pstar (and bfoldstar if
!present)
  function hbi_lnrhoreh_fromepsilon(w,eps1,eps2,eps3,Pstar,bfoldstar)     
    implicit none
    real(kp) :: hbi_lnrhoreh_fromepsilon
    real(kp), intent(in) :: w,eps1,eps2,eps3,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: n, mu
    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar
    real(kp) :: deltaNstar
    real(kp), dimension(3) :: hbistar

    logical, parameter :: printTest = .false.

    if (w.eq.1._kp/3._kp) then
       write(*,*)'hbi_lnrhoreh_fromepsilon: w = 1/3!'
       stop
    end if
    
    hbistar = hbi_xnmu_fromepsilon(eps1,eps2,eps3,bfoldstar)

    x = hbistar(1)
    n = hbistar(2)
    mu = hbistar(3)

    xEnd = hbi_x_endinf(n,mu)
    potEnd  = hbi_norm_potential(xEnd,n,mu)
    epsOneEnd = hbi_epsilon_one(xEnd,n,mu)
    
    potStar = hbi_norm_potential(x,n,mu)
    epsOneStar = hbi_epsilon_one(x,n,mu)

    if (.not.slowroll_validity(epsOneStar)) stop 'hbi_lnrhoreh_fromepsilon: slow-roll violated!'

   
    deltaNstar = hbi_efold_primitive(x,n,mu) - hbi_efold_primitive(xEnd,n,mu)

    if (printTest) then
       write(*,*)'eps1in= eps1comp= ',eps1, epsOneStar
       write(*,*)'eps2in= eps2comp= ',eps2,hbi_epsilon_two(x,n,mu)
       write(*,*)'eps3in= eps3comp= ',eps3,hbi_epsilon_three(x,n,mu)
    endif
    
    hbi_lnrhoreh_fromepsilon =  ln_rho_reheat(w,Pstar,epsOneStar,epsOneEnd,deltaNstar &
         ,potEnd/potStar)

  end function hbi_lnrhoreh_fromepsilon



end module hbireheat

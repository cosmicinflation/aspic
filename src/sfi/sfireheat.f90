!small field model reheating functions in the slow-roll approximations

module sfireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf,ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use sfisr, only : sfi_epsilon_one, sfi_epsilon_two, sfi_epsilon_three
  use sfisr, only : sfi_norm_potential
  use sfisr, only : sfi_x_endinf, sfi_efold_primitive
  implicit none

  private

  public sfi_x_star, sfi_lnrhoreh_max, sfi_lnrhoreh_fromepsilon, sfi_xpmu_fromepsilon
  public sfi_x_rrad, sfi_x_rreh


contains

!returns x given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the correspoding bfoldstar
  function sfi_x_star(p,mu,xend,w,lnRhoReh,Pstar,bfold)    
    implicit none
    real(kp) :: sfi_x_star
    real(kp), intent(in) :: p,mu,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: sfiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = sfi_epsilon_one(xEnd,p,mu)

    potEnd = sfi_norm_potential(xEnd,p,mu)
    primEnd = sfi_efold_primitive(xEnd,p,mu)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    sfiData%real1 = p
    sfiData%real2 = mu
    sfiData%real3 = w
    sfiData%real4 = calF + primEnd

    mini = 0._kp
    maxi = xEnd + epsilon(1._kp)

    x = zbrent(find_sfi_x_star,mini,maxi,tolFind,sfiData)
    sfi_x_star = x

    if (present(bfold)) then
       bfold = -(sfi_efold_primitive(x,p,mu) - primEnd)
    endif

    if (x.gt.1._kp) then
       if (display) write(*,*) 'sfi_x_star: phi>mu!'
    endif

  end function sfi_x_star

  function find_sfi_x_star(x,sfiData)   
    implicit none
    real(kp) :: find_sfi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: sfiData

    real(kp) :: primStar,p,mu,w,CalFplusPrimEnd,potStar,epsOneStar

    p=sfiData%real1
    mu = sfiData%real2
    w = sfiData%real3
    CalFplusPrimEnd = sfiData%real4

    primStar = sfi_efold_primitive(x,p,mu)
    epsOneStar = sfi_epsilon_one(x,p,mu)
    potStar = sfi_norm_potential(x,p,mu)


    find_sfi_x_star = find_reheat(PrimStar,calFplusPrimEnd,w,epsOneStar,potStar)

  end function find_sfi_x_star

!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function sfi_x_rrad(p,mu,xend,lnRrad,Pstar,bfold)    
    implicit none
    real(kp) :: sfi_x_rrad
    real(kp), intent(in) :: p,mu,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: sfiData
    
    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif
   
    epsOneEnd = sfi_epsilon_one(xEnd,p,mu)

    potEnd = sfi_norm_potential(xEnd,p,mu)
    primEnd = sfi_efold_primitive(xEnd,p,mu)
   
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    sfiData%real1 = p
    sfiData%real2 = mu
    sfiData%real3 = calF + primEnd

    mini = 0._kp
    maxi = xEnd + epsilon(1._kp)

    x = zbrent(find_sfi_x_rrad,mini,maxi,tolFind,sfiData)
    sfi_x_rrad = x

    if (present(bfold)) then
       bfold = -(sfi_efold_primitive(x,p,mu) - primEnd)
    endif

    if (x.gt.1._kp) then
       if (display) write(*,*) 'sfi_x_rrad: phi>mu!'
    endif

  end function sfi_x_rrad

  function find_sfi_x_rrad(x,sfiData)   
    implicit none
    real(kp) :: find_sfi_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: sfiData

    real(kp) :: primStar,p,mu,CalFplusPrimEnd,potStar,epsOneStar

    p=sfiData%real1
    mu = sfiData%real2
    CalFplusPrimEnd = sfiData%real3

    primStar = sfi_efold_primitive(x,p,mu)
    epsOneStar = sfi_epsilon_one(x,p,mu)
    potStar = sfi_norm_potential(x,p,mu)

    find_sfi_x_rrad = find_reheat_rrad(PrimStar,calFplusPrimEnd,epsOneStar,potStar)

  end function find_sfi_x_rrad



!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function sfi_x_rreh(p,mu,xend,lnRreh,bfold)    
    implicit none
    real(kp) :: sfi_x_rreh
    real(kp), intent(in) :: p,mu,xend,lnRreh
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: sfiData
    
    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif
   
    epsOneEnd = sfi_epsilon_one(xEnd,p,mu)

    potEnd = sfi_norm_potential(xEnd,p,mu)
    primEnd = sfi_efold_primitive(xEnd,p,mu)
   
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    sfiData%real1 = p
    sfiData%real2 = mu
    sfiData%real3 = calF + primEnd

    mini = 0._kp
    maxi = xEnd + epsilon(1._kp)

    x = zbrent(find_sfi_x_rreh,mini,maxi,tolFind,sfiData)
    sfi_x_rreh = x

    if (present(bfold)) then
       bfold = -(sfi_efold_primitive(x,p,mu) - primEnd)
    endif

    if (x.gt.1._kp) then
       if (display) write(*,*) 'sfi_x_rreh: phi>mu!'
    endif

  end function sfi_x_rreh

  function find_sfi_x_rreh(x,sfiData)   
    implicit none
    real(kp) :: find_sfi_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: sfiData

    real(kp) :: primStar,p,mu,CalFplusPrimEnd,potStar

    p=sfiData%real1
    mu = sfiData%real2
    CalFplusPrimEnd = sfiData%real3

    primStar = sfi_efold_primitive(x,p,mu)
    potStar = sfi_norm_potential(x,p,mu)

    find_sfi_x_rreh = find_reheat_rreh(PrimStar,calFplusPrimEnd,potStar)

  end function find_sfi_x_rreh



  function sfi_lnrhoreh_max(p,mu,xend,Pstar) 
    implicit none
    real(kp) :: sfi_lnrhoreh_max
    real(kp), intent(in) :: p,mu,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp), parameter :: wrad = 1._kp/3._kp
    real(kp), parameter :: junk= 0._kp
    real(kp) :: lnRhoEnd
    
    potEnd  = sfi_norm_potential(xEnd,p,mu)
    epsOneEnd = sfi_epsilon_one(xEnd,p,mu)
       
    x = sfi_x_star(p,mu,xend,wrad,junk,Pstar)    
    potStar = sfi_norm_potential(x,p,mu)
    epsOneStar = sfi_epsilon_one(x,p,mu)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'sfi_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    sfi_lnrhoreh_max = lnRhoEnd

  end function sfi_lnrhoreh_max

  
   
!return the unique x,p,mu giving eps123 (and bfold if input)
  function sfi_xpmu_fromepsilon(eps1,eps2,eps3,bfold)    
    implicit none
    real(kp), dimension(3) :: sfi_xpmu_fromepsilon
    real(kp), intent(in) :: eps1,eps2,eps3
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind = tolkp   
    type(transfert) :: sfiData
    real(kp) :: p,mu
    real(kp) :: x, xpowerp
    real(kp) :: xEnd

    if (.not.check_sfi_xpmu_fromepsilon(eps1,eps2,eps3)) then
       stop 'sfi_xpmu_fromepsilon: eps123 do not defined an unique small field model!'
    endif
    
    p = (2._kp/eps2)*( 8._kp*eps1**2 - 2._kp*eps1*eps2 + eps2**2 - eps2*eps3) &
         /(4._kp*eps1 + eps2 - 2._kp*eps3)
  
    xpowerp = 2._kp*eps1*(eps2 - 4._kp*eps1)/(eps2*(eps2 - eps3))

    x = xpowerp**(1._kp/p)

    mu = p/sqrt(2._kp*eps1) * x**(p-1)/(1._kp - xpowerp)
       

    if (present(bfold)) then        
       xEnd = sfi_x_endinf(p,mu)
       bfold = -(sfi_efold_primitive(x,p,mu)-sfi_efold_primitive(xEnd,p,mu))
    endif

    sfi_xpmu_fromepsilon(1) = x
    sfi_xpmu_fromepsilon(2) = p
    sfi_xpmu_fromepsilon(3) = mu

  end function sfi_xpmu_fromepsilon

!check if eps123 uniquely define x,mu,p
  function check_sfi_xpmu_fromepsilon(eps1,eps2,eps3)
    real(kp), intent(in) :: eps1,eps2,eps3
    logical :: check_sfi_xpmu_fromepsilon

    real(kp) :: xpowerp, p

    check_sfi_xpmu_fromepsilon = .false.
    
    if (eps2.lt.0._kp) return

    if ((4._kp*eps1 + eps2 - 2._kp*eps3).eq.0._kp) then      
       if (display) write(*,*)'check_sfi_xpmu_fromepsilon: p -> infty!'
       return
    endif

    xpowerp = 2._kp*eps1*(eps2 - 4._kp*eps1)/(eps2*(eps2 - eps3))

    if ((xpowerp.gt.1._kp).or.(xpowerp.lt.0._kp)) then
       if (display) write(*,*)'check_sfi_xpmu_frompepsilon: x^p= ',xpowerp
       return       
    endif

    p = (2._kp/eps2)*( 8._kp*eps1**2 - 2._kp*eps1*eps2 + eps2**2 - eps2*eps3) &
         /(4._kp*eps1 + eps2 - 2._kp*eps3)

    if (p.lt.1._kp) then
       if (display) write(*,*)'check_sfi_xpmu_fromepsilon: p<1'
       return
    endif

    check_sfi_xpmu_fromepsilon = .true.

  end function check_sfi_xpmu_fromepsilon


!returns lnrhoreh from eps123, wreh and Pstar (and bfoldstar if
!present)
  function sfi_lnrhoreh_fromepsilon(w,eps1,eps2,eps3,Pstar,bfoldstar)     
    implicit none
    real(kp) :: sfi_lnrhoreh_fromepsilon
    real(kp), intent(in) :: w,eps1,eps2,eps3,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: p, mu
    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar
    real(kp) :: deltaNstar
    real(kp), dimension(3) :: sfistar

    logical, parameter :: printTest = .false.

    if (w.eq.1._kp/3._kp) then
       write(*,*)'sfi_lnrhoreh_fromepsilon: w = 1/3!'
       stop
    end if
    
    sfistar = sfi_xpmu_fromepsilon(eps1,eps2,eps3,bfoldstar)

    x = sfistar(1)
    p = sfistar(2)
    mu = sfistar(3)

    xEnd = sfi_x_endinf(p,mu)       
    potEnd  = sfi_norm_potential(xEnd,p,mu)
    epsOneEnd = sfi_epsilon_one(xEnd,p,mu)
    
    potStar = sfi_norm_potential(x,p,mu)
    epsOneStar = sfi_epsilon_one(x,p,mu)

    if (.not.slowroll_validity(epsOneStar)) stop 'sfi_lnrhoreh_fromepsilon: slow-roll violated!'

   
    deltaNstar = sfi_efold_primitive(x,p,mu) - sfi_efold_primitive(xEnd,p,mu)

    if (printTest) then
       write(*,*)'eps1in= eps1comp= ',eps1, epsOneStar
       write(*,*)'eps2in= eps2comp= ',eps2,sfi_epsilon_two(x,p,mu)
       write(*,*)'eps3in= eps3comp= ',eps3,sfi_epsilon_three(x,p,mu)
    endif
    
    sfi_lnrhoreh_fromepsilon =  ln_rho_reheat(w,Pstar,epsOneStar,epsOneEnd,deltaNstar &
         ,potEnd/potStar)

  end function sfi_lnrhoreh_fromepsilon



end module sfireheat

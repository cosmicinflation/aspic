!large field model reheating functions in the slow-roll approximations

module lfireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : get_calfconst_rrad, find_reheat_rrad
  use srreheat, only : get_calfconst_rreh, find_reheat_rreh
  use lfisr, only : lfi_epsilon_one, lfi_epsilon_two, lfi_epsilon_three
  use lfisr, only : lfi_norm_potential
  use lfisr, only : lfi_x_endinf, lfi_efold_primitive
  implicit none

  private

  public lfi_x_star, lfi_lnrhoreh_max
  public lfi_xp_fromepsilon, lfi_lnrhoreh_fromepsilon 

  public lfi_x_rrad, lfi_x_rreh

contains

!returns x given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function lfi_x_star(p,xend,w,lnRhoReh,Pstar,bfoldstar)   
    implicit none
    real(kp) :: lfi_x_star
    real(kp), intent(in) :: p,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: lfiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = lfi_epsilon_one(xEnd,p)
    potEnd = lfi_norm_potential(xEnd,p)
    primEnd = lfi_efold_primitive(xEnd,p)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    lfiData%real1 = p    
    lfiData%real2 = w
    lfiData%real3 = calF + primEnd

    mini = xEnd
    maxi = 1._kp/epsilon(1._kp)

    x = zbrent(find_lfi_x_star,mini,maxi,tolzbrent,lfiData)
    lfi_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (lfi_efold_primitive(x,p) - primEnd)
    endif

  end function lfi_x_star

  function find_lfi_x_star(x,lfiData)   
    implicit none
    real(kp) :: find_lfi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: lfiData

    real(kp) :: primStar,p,w,CalFplusprimEnd,potStar
    real(kp) :: epsOneStar

    p=lfiData%real1
    w = lfiData%real2
    CalFplusprimEnd = lfiData%real3

    primStar = lfi_efold_primitive(x,p)
    epsOneStar = lfi_epsilon_one(x,p)    
    potStar = lfi_norm_potential(x,p)

    find_lfi_x_star = find_reheat(primStar,calFplusprimEnd,w&
         ,epsOneStar,potStar)
  
  end function find_lfi_x_star

!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function lfi_x_rrad(p,xend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: lfi_x_rrad
    real(kp), intent(in) :: p,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: lfiData
    

    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for instantaneous reheating!'
    endif
    
    epsOneEnd = lfi_epsilon_one(xEnd,p)
    potEnd = lfi_norm_potential(xEnd,p)
    primEnd = lfi_efold_primitive(xEnd,p)
   
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    lfiData%real1 = p    
    lfiData%real2 = calF + primEnd

    mini = xEnd
    maxi = 1._kp/epsilon(1._kp)

    x = zbrent(find_lfi_x_rrad,mini,maxi,tolzbrent,lfiData)
    lfi_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (lfi_efold_primitive(x,p) - primEnd)
    endif

  end function lfi_x_rrad

  function find_lfi_x_rrad(x,lfiData)   
    implicit none
    real(kp) :: find_lfi_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: lfiData

    real(kp) :: primStar,p,w,CalFplusprimEnd
    real(kp) :: potStar
    real(kp) :: epsOneStar

    p=lfiData%real1    
    CalFplusprimEnd = lfiData%real2

    primStar = lfi_efold_primitive(x,p)
    epsOneStar = lfi_epsilon_one(x,p)
    potStar = lfi_norm_potential(x,p)

    find_lfi_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd &
         ,epsOneStar,potStar)
  
  end function find_lfi_x_rrad


!returns x given potential parameters and lnR (no need of Pstar, lnR
!is optimal for CMB). If present, returns the corresponding bfoldstar
  function lfi_x_rreh(p,xend,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: lfi_x_rreh
    real(kp), intent(in) :: p,xend,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: lfiData
    

    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for instantaneous reheating!'
    endif
    
    epsOneEnd = lfi_epsilon_one(xEnd,p)
    potEnd = lfi_norm_potential(xEnd,p)
    primEnd = lfi_efold_primitive(xEnd,p)
   
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    lfiData%real1 = p    
    lfiData%real2 = calF + primEnd

    mini = xEnd
    maxi = 1._kp/epsilon(1._kp)

    x = zbrent(find_lfi_x_rreh,mini,maxi,tolzbrent,lfiData)
    lfi_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (lfi_efold_primitive(x,p) - primEnd)
    endif

  end function lfi_x_rreh

  function find_lfi_x_rreh(x,lfiData)   
    implicit none
    real(kp) :: find_lfi_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: lfiData

    real(kp) :: primStar,p,w,CalFplusprimEnd
    real(kp) :: potStar

    p=lfiData%real1    
    CalFplusprimEnd = lfiData%real2

    primStar = lfi_efold_primitive(x,p)
    potStar = lfi_norm_potential(x,p)

    find_lfi_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd &
         ,potStar)
  
  end function find_lfi_x_rreh




  function lfi_lnrhoreh_max(p,xend,Pstar) 
    implicit none
    real(kp) :: lfi_lnrhoreh_max
    real(kp), intent(in) :: p,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd

    
    potEnd  = lfi_norm_potential(xEnd,p)
    epsOneEnd = lfi_epsilon_one(xEnd,p)

!   Trick to return x such that rho_reh=rho_end

    x = lfi_x_star(p,xend,wrad,junk,Pstar)    
    potStar = lfi_norm_potential(x,p)
    epsOneStar = lfi_epsilon_one(x,p)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'lfi_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    lfi_lnrhoreh_max = lnRhoEnd


  end function lfi_lnrhoreh_max

  
   
!return the unique p,x giving eps12 (and bfold if input)
  function lfi_xp_fromepsilon(eps1,eps2,bfold)    
    implicit none
    real(kp), dimension(2) :: lfi_xp_fromepsilon
    real(kp), intent(in) :: eps1,eps2
    real(kp), intent(out), optional :: bfold
    
    real(kp) :: x, xEnd, p

    if (eps2.le.0._kp) then
       stop 'lfi_xp_fromepsilon: eps2<=0'
    endif

    p = 4._kp*eps1/eps2
    x = sqrt(8._kp*eps1/eps2/eps2)

    lfi_xp_fromepsilon(1) = x
    lfi_xp_fromepsilon(2) = p
    
    xEnd = lfi_x_endinf(p)    
    
    if (present(bfold)) then        
       bfold = -(lfi_efold_primitive(x,p)-lfi_efold_primitive(xEnd,p))
    endif

  end function lfi_xp_fromepsilon

  
!returns lnrhoreh from eps12, wreh and Pstar (and bfoldstar if
!present)
  function lfi_lnrhoreh_fromepsilon(w,eps1,eps2,Pstar,bfoldstar)     
    implicit none
    real(kp) :: lfi_lnrhoreh_fromepsilon
    real(kp), intent(inout) :: w
    real(kp), intent(in) :: eps1,eps2,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: p
    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar
    real(kp) :: deltaNstar
    real(kp), dimension(2) :: lfistar

    logical, parameter :: printTest = .false.
    logical, parameter :: enforceWofp = .true.

    lfistar = lfi_xp_fromepsilon(eps1,eps2,bfoldstar)    
    x = lfistar(1)
    p = lfistar(2)

    if (enforceWofp) then
       w = (p-2._kp)/(p+2._kp)
    endif

    xEnd = lfi_x_endinf(p)       
    potEnd  = lfi_norm_potential(xEnd,p)
    epsOneEnd = lfi_epsilon_one(xEnd,p)
    
    potStar = lfi_norm_potential(x,p)
    epsOneStar = lfi_epsilon_one(x,p)

    if (.not.slowroll_validity(epsOneStar)) stop 'lfi_lnrhoreh: cannot trust slow-roll!'
                   
    deltaNstar = lfi_efold_primitive(x,p) - lfi_efold_primitive(xEnd,p)

    if (printTest) then
       write(*,*)'eps1in= eps1comp= ',eps1, epsOneStar
       write(*,*)'eps2in= eps2comp= ',eps2,lfi_epsilon_two(x,p)
       write(*,*)'eps3in= eps3comp= ',eps2,lfi_epsilon_three(x,p)
    endif
           
    lfi_lnrhoreh_fromepsilon =  ln_rho_reheat(w,Pstar,epsOneStar,epsOneEnd,deltaNstar &
         ,potEnd/potStar)
    
  end function lfi_lnrhoreh_fromepsilon



end module lfireheat

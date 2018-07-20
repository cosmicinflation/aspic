!radion gauge model reheating functions in the slow-roll approximations

module rgireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use rgisr, only : rgi_epsilon_one, rgi_epsilon_two, rgi_epsilon_three
  use rgisr, only : rgi_norm_potential, rgi_numacc_x_epsonenull
  use rgisr, only : rgi_x_endinf, rgi_efold_primitive
  implicit none

  private

  public rgi_x_star, rgi_lnrhoreh_max 
  public rgi_xp_fromepsilon, rgi_lnrhoreh_fromepsilon 
  public rgi_x_rrad, rgi_x_rreh

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function rgi_x_star(alpha,xend,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: rgi_x_star
    real(kp), intent(in) :: alpha,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: rgiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = rgi_epsilon_one(xEnd,alpha)
    potEnd = rgi_norm_potential(xEnd,alpha)
    primEnd = rgi_efold_primitive(xEnd,alpha)

   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    rgiData%real1 = alpha    
    rgiData%real2 = w
    rgiData%real3 = calF + primEnd

    mini = xEnd
    maxi = rgi_numacc_x_epsonenull(alpha)

    x = zbrent(find_rgi_x_star,mini,maxi,tolzbrent,rgiData)
    rgi_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (rgi_efold_primitive(x,alpha) - primEnd)
    endif

  end function rgi_x_star

  function find_rgi_x_star(x,rgiData)   
    implicit none
    real(kp) :: find_rgi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rgiData

    real(kp) :: primStar,alpha,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=rgiData%real1
    w = rgiData%real2
    CalFplusprimEnd = rgiData%real3

    primStar = rgi_efold_primitive(x,alpha)
    epsOneStar = rgi_epsilon_one(x,alpha)
    potStar = rgi_norm_potential(x,alpha)

    find_rgi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_rgi_x_star

!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function rgi_x_rrad(alpha,xend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: rgi_x_rrad
    real(kp), intent(in) :: alpha,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: rgiData
    

    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif
        
    epsOneEnd = rgi_epsilon_one(xEnd,alpha)
    potEnd = rgi_norm_potential(xEnd,alpha)
    primEnd = rgi_efold_primitive(xEnd,alpha)

   
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    rgiData%real1 = alpha    
    rgiData%real2 = calF + primEnd

    mini = xEnd
    maxi = rgi_numacc_x_epsonenull(alpha)

    x = zbrent(find_rgi_x_rrad,mini,maxi,tolzbrent,rgiData)
    rgi_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (rgi_efold_primitive(x,alpha) - primEnd)
    endif

  end function rgi_x_rrad

  function find_rgi_x_rrad(x,rgiData)   
    implicit none
    real(kp) :: find_rgi_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rgiData

    real(kp) :: primStar,alpha,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=rgiData%real1
    CalFplusprimEnd = rgiData%real2

    primStar = rgi_efold_primitive(x,alpha)
    epsOneStar = rgi_epsilon_one(x,alpha)
    potStar = rgi_norm_potential(x,alpha)

    find_rgi_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)
  
  end function find_rgi_x_rrad



!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function rgi_x_rreh(alpha,xend,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: rgi_x_rreh
    real(kp), intent(in) :: alpha,xend,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: rgiData
    

    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif
        
    epsOneEnd = rgi_epsilon_one(xEnd,alpha)
    potEnd = rgi_norm_potential(xEnd,alpha)
    primEnd = rgi_efold_primitive(xEnd,alpha)

   
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    rgiData%real1 = alpha    
    rgiData%real2 = calF + primEnd

    mini = xEnd
    maxi = rgi_numacc_x_epsonenull(alpha)

    x = zbrent(find_rgi_x_rreh,mini,maxi,tolzbrent,rgiData)
    rgi_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (rgi_efold_primitive(x,alpha) - primEnd)
    endif

  end function rgi_x_rreh

  function find_rgi_x_rreh(x,rgiData)   
    implicit none
    real(kp) :: find_rgi_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rgiData

    real(kp) :: primStar,alpha,w,CalFplusprimEnd,potStar

    alpha=rgiData%real1
    CalFplusprimEnd = rgiData%real2

    primStar = rgi_efold_primitive(x,alpha)
    potStar = rgi_norm_potential(x,alpha)

    find_rgi_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)
  
  end function find_rgi_x_rreh



  function rgi_lnrhoreh_max(alpha,xend,Pstar) 
    implicit none
    real(kp) :: rgi_lnrhoreh_max
    real(kp), intent(in) :: alpha,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd
    
    potEnd  = rgi_norm_potential(xEnd,alpha)
    epsOneEnd = rgi_epsilon_one(xEnd,alpha)


!   Trick to return x such that rho_reh=rho_end

    x = rgi_x_star(alpha,xend,wrad,junk,Pstar)    
    potStar = rgi_norm_potential(x,alpha)
    epsOneStar = rgi_epsilon_one(x,alpha)

    
    if (.not.slowroll_validity(epsOneStar)) stop 'rgi_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    rgi_lnrhoreh_max = lnRhoEnd

  end function rgi_lnrhoreh_max

  
   
!return the unique alpha,x giving eps12 (and bfold if input)
  function rgi_xp_fromepsilon(eps1,eps2,bfold)    
    implicit none
    real(kp), dimension(2) :: rgi_xp_fromepsilon
    real(kp), intent(in) :: eps1,eps2
    real(kp), intent(out), optional :: bfold
    
    real(kp) :: x, xEnd, alpha

    if (eps2.le.0._kp) then
       stop 'rgi_xp_fromepsilon: eps2<=0'
    endif

    if (eps1.le.0._kp) then
       stop 'rgi_xp_fromepsilon: eps1<=0'
    endif

    alpha = 432._kp*eps1**2/((eps2-2._kp*eps1)*(4._kp*eps1+eps2)**2)
    x = 6._kp*sqrt(2._kp*eps1)/(4._kp*eps1+eps2)

    rgi_xp_fromepsilon(1) = x
    rgi_xp_fromepsilon(2) = alpha
    
    xEnd = rgi_x_endinf(alpha)    
    
    if (present(bfold)) then        
       bfold = -(rgi_efold_primitive(x,alpha)-rgi_efold_primitive(xEnd,alpha))
    endif

  end function rgi_xp_fromepsilon

  
!returns lnrhoreh from eps12, wreh and Pstar (and bfoldstar if
!present)
  function rgi_lnrhoreh_fromepsilon(w,eps1,eps2,Pstar,bfoldstar)     
    implicit none
    real(kp) :: rgi_lnrhoreh_fromepsilon
    real(kp), intent(inout) :: w
    real(kp), intent(in) :: eps1,eps2,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: alpha
    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar
    real(kp) :: deltaNstar
    real(kp), dimension(2) :: rgistar

    logical, parameter :: printTest = .false.

    rgistar = rgi_xp_fromepsilon(eps1,eps2,bfoldstar)    
    x = rgistar(1)
    alpha = rgistar(2)

    w=0._kp

    xEnd = rgi_x_endinf(alpha)       
    potEnd  = rgi_norm_potential(xEnd,alpha)
    epsOneEnd = rgi_epsilon_one(xEnd,alpha)
    
    potStar = rgi_norm_potential(x,alpha)
    epsOneStar = rgi_epsilon_one(x,alpha)

    if (.not.slowroll_validity(epsOneStar)) stop 'rgi_lnrhoreh: cannot trust slow-roll!'
                   
    deltaNstar = rgi_efold_primitive(x,alpha) - rgi_efold_primitive(xEnd,alpha)

    if (printTest) then
       write(*,*)'eps1in= eps1comp= ',eps1, epsOneStar
       write(*,*)'eps2in= eps2comp= ',eps2,rgi_epsilon_two(x,alpha)
       write(*,*)'eps3in= eps3comp= ',eps2,rgi_epsilon_three(x,alpha)
    endif
           
    rgi_lnrhoreh_fromepsilon =  ln_rho_reheat(w,Pstar,epsOneStar,epsOneEnd,deltaNstar &
         ,potEnd/potStar)
    
  end function rgi_lnrhoreh_fromepsilon



end module rgireheat

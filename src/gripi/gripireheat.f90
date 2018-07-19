!GRIP inflation reheating functions in the slow-roll approximations

module gripireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use gripisr, only : gripi_epsilon_one, gripi_epsilon_two, gripi_epsilon_three
  use gripisr, only : gripi_norm_potential, gripi_x_endinf, gripi_x_epsonemin
  use gripisr, only : gripi_efold_primitive
  implicit none

  private

  public gripi_x_star, gripi_lnrhoreh_max 
  public gripi_x_rrad, gripi_x_rreh

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function gripi_x_star(alpha,phi0,xend,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: gripi_x_star
    real(kp), intent(in) :: alpha,phi0,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: gripiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = gripi_epsilon_one(xEnd,alpha,phi0)
    potEnd = gripi_norm_potential(xEnd,alpha,phi0)
    primEnd = gripi_efold_primitive(xEnd,alpha,phi0)

    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    gripiData%real1 = alpha 
    gripiData%real2 = phi0 
    gripiData%real3 = w
    gripiData%real4 = calF + primEnd

    mini = xend
    maxi = gripi_x_epsonemin(alpha,phi0) - epsilon(1._kp)

    x = zbrent(find_gripi_x_star,mini,maxi,tolzbrent,gripiData)
    gripi_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (gripi_efold_primitive(x,alpha,phi0) - primEnd)
    endif


  end function gripi_x_star

  function find_gripi_x_star(x,gripiData)   
    implicit none
    real(kp) :: find_gripi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: gripiData

    real(kp) :: primStar,alpha,phi0,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=gripiData%real1
    phi0=gripiData%real2
    w = gripiData%real3
    CalFplusprimEnd = gripiData%real4

    primStar = gripi_efold_primitive(x,alpha,phi0)
    epsOneStar = gripi_epsilon_one(x,alpha,phi0)
    potStar = gripi_norm_potential(x,alpha,phi0)

    find_gripi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)

  end function find_gripi_x_star

!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function gripi_x_rrad(alpha,phi0,xend,lnRRad,Pstar,bfoldstar)        
    implicit none
    real(kp) :: gripi_x_rrad
    real(kp), intent(in) :: alpha,phi0,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: gripiData
    
    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = gripi_epsilon_one(xEnd,alpha,phi0)
    potEnd = gripi_norm_potential(xEnd,alpha,phi0)
    primEnd = gripi_efold_primitive(xEnd,alpha,phi0)

    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)


    gripiData%real1 = alpha
    gripiData%real2 = phi0
    gripiData%real3 = calF + primEnd

    mini = xend
    maxi = gripi_x_epsonemin(alpha,phi0) - epsilon(1._kp)
    
    x = zbrent(find_gripi_x_rrad,mini,maxi,tolzbrent,gripiData)
    gripi_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (gripi_efold_primitive(x,alpha,phi0) - primEnd)
    endif

  end function gripi_x_rrad

  function find_gripi_x_rrad(x,gripiData)   
    implicit none
    real(kp) :: find_gripi_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: gripiData

    real(kp) :: primStar,alpha,phi0,CalFplusprimEnd,potStar,epsOneStar

    alpha = gripiData%real1
    phi0 = gripiData%real2
    CalFplusprimEnd = gripiData%real3

    primStar = gripi_efold_primitive(x,alpha,phi0)
    epsOneStar = gripi_epsilon_one(x,alpha,phi0)
    potStar = gripi_norm_potential(x,alpha,phi0)

    find_gripi_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd &
         ,epsOneStar,potStar)
  
  end function find_gripi_x_rrad


  !returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function gripi_x_rreh(alpha,phi0,xend,lnRReh,bfoldstar)    
    implicit none
    real(kp) :: gripi_x_rreh
    real(kp), intent(in) :: alpha,phi0,xend,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: gripiData
    

    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = gripi_epsilon_one(xEnd,alpha,phi0)
    potEnd = gripi_norm_potential(xEnd,alpha,phi0)
    primEnd = gripi_efold_primitive(xEnd,alpha,phi0)

    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    gripiData%real1 = alpha
    gripiData%real2 = phi0
    gripiData%real3 = calF + primEnd

    mini = xend
    maxi = gripi_x_epsonemin(alpha,phi0) - epsilon(1._kp)
    
    x = zbrent(find_gripi_x_rreh,mini,maxi,tolzbrent,gripiData)
    gripi_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (gripi_efold_primitive(x,alpha,phi0) - primEnd)
    endif

  end function gripi_x_rreh

  function find_gripi_x_rreh(x,gripiData)   
    implicit none
    real(kp) :: find_gripi_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: gripiData

    real(kp) :: primStar,alpha,phi0,CalFplusprimEnd,potStar

    alpha = gripiData%real1
    phi0 = gripiData%real2
    CalFplusprimEnd = gripiData%real3

    primStar = gripi_efold_primitive(x,alpha,phi0)    
    potStar = gripi_norm_potential(x,alpha,phi0)

    find_gripi_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd &
         ,potStar)
  
  end function find_gripi_x_rreh



  function gripi_lnrhoreh_max(alpha,phi0,xend,Pstar) 
    implicit none
    real(kp) :: gripi_lnrhoreh_max
    real(kp), intent(in) :: alpha,phi0,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd
    
    potEnd  = gripi_norm_potential(xEnd,alpha,phi0)
    epsOneEnd = gripi_epsilon_one(xEnd,alpha,phi0)

!   Trick to return x such that rho_reh=rho_end

    x = gripi_x_star(alpha,phi0,xend,wrad,junk,Pstar)    
    potStar = gripi_norm_potential(x,alpha,phi0)
    epsOneStar = gripi_epsilon_one(x,alpha,phi0)

    
!    if (.not.slowroll_validity(epsOneStar)) stop 'gripi_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    gripi_lnrhoreh_max = lnRhoEnd

  end function gripi_lnrhoreh_max

  
end module gripireheat

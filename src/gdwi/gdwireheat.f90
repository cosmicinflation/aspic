!double well inflation reheating functions in the slow-roll approximations

module gdwireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use gdwisr, only : gdwi_epsilon_one, gdwi_epsilon_two, gdwi_epsilon_three
  use gdwisr, only : gdwi_norm_potential
  use gdwisr, only : gdwi_x_endinf, gdwi_efold_primitive
  implicit none

  private

  public gdwi_x_star, gdwi_lnrhoreh_max
  public gdwi_x_rrad, gdwi_x_rreh
  

contains

!returns x =phi/phi0 such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function gdwi_x_star(p,phi0,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: gdwi_x_star
    real(kp), intent(in) :: p,phi0,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: gdwiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = gdwi_x_endinf(p,phi0)

    epsOneEnd = gdwi_epsilon_one(xEnd,p,phi0)
    potEnd = gdwi_norm_potential(xEnd,p,phi0)

    primEnd = gdwi_efold_primitive(xEnd,p,phi0)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    gdwiData%real1 = p
    gdwiData%real2 = phi0
    gdwiData%real3 = w
    gdwiData%real4 = calF + primEnd

    mini = epsilon(1._kp)

    maxi = xEnd

   
    x = zbrent(find_gdwi_x_star,mini,maxi,tolzbrent,gdwiData)
    gdwi_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (gdwi_efold_primitive(x,p,phi0) - primEnd)
    endif

  end function gdwi_x_star


  function find_gdwi_x_star(x,gdwiData)   
    implicit none
    real(kp) :: find_gdwi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: gdwiData

    real(kp) :: primStar,p,phi0,w,CalFplusprimEnd,potStar,epsOneStar

    p = gdwiData%real1
    phi0=gdwiData%real2
    w = gdwiData%real3
    CalFplusprimEnd = gdwiData%real4

    primStar = gdwi_efold_primitive(x,p,phi0)
    epsOneStar = gdwi_epsilon_one(x,p,phi0)
    potStar = gdwi_norm_potential(x,p,phi0)

    find_gdwi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_gdwi_x_star



!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function gdwi_x_rrad(p,phi0,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: gdwi_x_rrad
    real(kp), intent(in) :: p,phi0,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: gdwiData
    

    if (lnRRad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif

    xEnd = gdwi_x_endinf(p,phi0)

    epsOneEnd = gdwi_epsilon_one(xEnd,p,phi0)
    potEnd = gdwi_norm_potential(xEnd,p,phi0)

    primEnd = gdwi_efold_primitive(xEnd,p,phi0)
   
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    gdwiData%real1 = p
    gdwiData%real2 = phi0
    gdwiData%real3 = calF + primEnd

    mini = epsilon(1._kp)
    maxi = xEnd
   
    x = zbrent(find_gdwi_x_rrad,mini,maxi,tolzbrent,gdwiData)
    gdwi_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (gdwi_efold_primitive(x,p,phi0) - primEnd)
    endif

  end function gdwi_x_rrad


  function find_gdwi_x_rrad(x,gdwiData)   
    implicit none
    real(kp) :: find_gdwi_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: gdwiData

    real(kp) :: primStar,p,phi0,CalFplusprimEnd,potStar,epsOneStar

    p = gdwiData%real1
    phi0=gdwiData%real2
    CalFplusprimEnd = gdwiData%real3

    primStar = gdwi_efold_primitive(x,p,phi0)
    epsOneStar = gdwi_epsilon_one(x,p,phi0)
    potStar = gdwi_norm_potential(x,p,phi0)

    find_gdwi_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)
  
  end function find_gdwi_x_rrad


  
!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function gdwi_x_rreh(p,phi0,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: gdwi_x_rreh
    real(kp), intent(in) :: p,phi0,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: gdwiData
    

    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif

    xEnd = gdwi_x_endinf(p,phi0)

    epsOneEnd = gdwi_epsilon_one(xEnd,p,phi0)
    potEnd = gdwi_norm_potential(xEnd,p,phi0)

    primEnd = gdwi_efold_primitive(xEnd,p,phi0)
   
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    gdwiData%real1 = p
    gdwiData%real2 = phi0
    gdwiData%real3 = calF + primEnd

    mini = epsilon(1._kp)
    maxi = xEnd
   
    x = zbrent(find_gdwi_x_rreh,mini,maxi,tolzbrent,gdwiData)
    gdwi_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (gdwi_efold_primitive(x,p,phi0) - primEnd)
    endif

  end function gdwi_x_rreh


  function find_gdwi_x_rreh(x,gdwiData)   
    implicit none
    real(kp) :: find_gdwi_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: gdwiData

    real(kp) :: primStar,p,phi0,CalFplusprimEnd,potStar

    p = gdwiData%real1
    phi0=gdwiData%real2
    CalFplusprimEnd = gdwiData%real3

    primStar = gdwi_efold_primitive(x,p,phi0)    
    potStar = gdwi_norm_potential(x,p,phi0)

    find_gdwi_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)
  
  end function find_gdwi_x_rreh



  function gdwi_lnrhoreh_max(p,phi0,Pstar) 
    implicit none
    real(kp) :: gdwi_lnrhoreh_max
    real(kp), intent(in) :: p,phi0,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd
    
    xEnd = gdwi_x_endinf(p,phi0)


    potEnd  = gdwi_norm_potential(xEnd,p,phi0)

    epsOneEnd = gdwi_epsilon_one(xEnd,p,phi0)


!   Trick to return x such that rho_reh=rho_end

    x = gdwi_x_star(p,phi0,wrad,junk,Pstar)  

 
    potStar = gdwi_norm_potential(x,p,phi0)
    epsOneStar = gdwi_epsilon_one(x,p,phi0)

   ! PRINT*,'gdwi_lnrhoreh_max   :xstar=',x,'  potStar=',potStar,'  epsOneStar=',epsOneStar

    
    if (.not.slowroll_validity(epsOneStar)) stop 'gdwi_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    gdwi_lnrhoreh_max = lnRhoEnd

  end function gdwi_lnrhoreh_max

  
end module gdwireheat

!MSSM inflation reheating functions in the slow-roll approximations

module gmssmireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use gmssmisr, only : gmssmi_epsilon_one, gmssmi_epsilon_two, gmssmi_epsilon_three
  use gmssmisr, only : gmssmi_norm_potential, gmssmi_x_endinf, gmssmi_x_epsonemin
  use gmssmisr, only : gmssmi_efold_primitive
  implicit none

  private

  public gmssmi_x_star, gmssmi_lnrhoreh_max 

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function gmssmi_x_star(alpha,phi0,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: gmssmi_x_star
    real(kp), intent(in) :: alpha,phi0,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: gmssmiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = gmssmi_x_endinf(alpha,phi0)
    epsOneEnd = gmssmi_epsilon_one(xEnd,alpha,phi0)
    potEnd = gmssmi_norm_potential(xEnd,alpha,phi0)
    primEnd = gmssmi_efold_primitive(xEnd,alpha,phi0)

    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    gmssmiData%real1 = alpha 
    gmssmiData%real2 = phi0 
    gmssmiData%real3 = xEnd
    gmssmiData%real4 = w
    gmssmiData%real5 = calF + primEnd

    mini = xend

    if (alpha .lt. 1._kp) then
	maxi = gmssmi_x_epsonemin(alpha)*(1._kp-100000._kp*epsilon(1._kp))
    else
	maxi = gmssmi_x_epsonemin(alpha)*(1._kp-100000._kp*epsilon(1._kp)) !local maximum of the potential
    endif

    x = zbrent(find_gmssmi_x_star,mini,maxi,tolzbrent,gmssmiData)
    gmssmi_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (gmssmi_efold_primitive(x,alpha,phi0) - primEnd)
    endif


  end function gmssmi_x_star

  function find_gmssmi_x_star(x,gmssmiData)   
    implicit none
    real(kp) :: find_gmssmi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: gmssmiData

    real(kp) :: primStar,alpha,phi0,xEnd,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=gmssmiData%real1
    phi0=gmssmiData%real2
    xEnd=gmssmiData%real3
    w = gmssmiData%real4
    CalFplusprimEnd = gmssmiData%real5

    primStar = gmssmi_efold_primitive(x,alpha,phi0)
    epsOneStar = gmssmi_epsilon_one(x,alpha,phi0)
    potStar = gmssmi_norm_potential(x,alpha,phi0)

    find_gmssmi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)

  
  end function find_gmssmi_x_star



  function gmssmi_lnrhoreh_max(alpha,phi0,Pstar) 
    implicit none
    real(kp) :: gmssmi_lnrhoreh_max
    real(kp), intent(in) :: alpha,phi0,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = gmssmi_x_endinf(alpha,phi0)
    potEnd  = gmssmi_norm_potential(xEnd,alpha,phi0)
    epsOneEnd = gmssmi_epsilon_one(xEnd,alpha,phi0)

!   Trick to return x such that rho_reh=rho_end

    x = gmssmi_x_star(alpha,phi0,wrad,junk,Pstar)    
    potStar = gmssmi_norm_potential(x,alpha,phi0)
    epsOneStar = gmssmi_epsilon_one(x,alpha,phi0)

    
!    if (.not.slowroll_validity(epsOneStar)) stop 'gmssmi_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    gmssmi_lnrhoreh_max = lnRhoEnd

  end function gmssmi_lnrhoreh_max

  
end module gmssmireheat

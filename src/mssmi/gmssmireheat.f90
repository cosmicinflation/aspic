!MSSM inflation reheating functions in the slow-roll approximations

module gmssmireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use gmssmisr, only : gmssmi_epsilon_one, gmssmi_epsilon_two, gmssmi_epsilon_three
  use gmssmisr, only : gmssmi_norm_potential, gmssmi_x_endinf, gmssmi_x_epsonemin
  use gmssmisr, only :  gmssmi_efold_primitive
  implicit none

  private

  public gmssmi_x_star, gmssmi_lnrhoend 

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function gmssmi_x_star(alpha,beta,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: gmssmi_x_star
    real(kp), intent(in) :: alpha,beta,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: gmssmiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = gmssmi_x_endinf(alpha,beta)
    epsOneEnd = gmssmi_epsilon_one(xEnd,alpha,beta)
    potEnd = gmssmi_norm_potential(xEnd,alpha,beta)
    primEnd = gmssmi_efold_primitive(xEnd,alpha,beta)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    gmssmiData%real1 = alpha 
    gmssmiData%real2 = beta 
    gmssmiData%real3 = xEnd
    gmssmiData%real4 = w
    gmssmiData%real5 = calF + primEnd

    mini = xEnd

    if (alpha**2/beta>20._kp/9._kp) then
	maxi = gmssmi_x_epsonemin(alpha,beta)!*(1._kp-1.*epsilon(1._kp)) !local maximum of the potential
    else
	maxi=100._kp*gmssmi_x_epsonemin(alpha,beta)
    endif

    x = zbrent(find_gmssmi_x_star,mini,maxi,tolzbrent,gmssmiData)
    gmssmi_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (gmssmi_efold_primitive(x,alpha,beta) - primEnd)
    endif

  end function gmssmi_x_star

  function find_gmssmi_x_star(x,gmssmiData)   
    implicit none
    real(kp) :: find_gmssmi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: gmssmiData

    real(kp) :: primStar,alpha,beta,xEnd,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=gmssmiData%real1
    beta=gmssmiData%real2
    xEnd=gmssmiData%real3
    w = gmssmiData%real4
    CalFplusprimEnd = gmssmiData%real5

    primStar = gmssmi_efold_primitive(x,alpha,beta)
    epsOneStar = gmssmi_epsilon_one(x,alpha,beta)
    potStar = gmssmi_norm_potential(x,alpha,beta)

    find_gmssmi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_gmssmi_x_star



  function gmssmi_lnrhoend(alpha,beta,Pstar) 
    implicit none
    real(kp) :: gmssmi_lnrhoend
    real(kp), intent(in) :: alpha,beta,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = gmssmi_x_endinf(alpha,beta)
    potEnd  = gmssmi_norm_potential(xEnd,alpha,beta)
    epsOneEnd = gmssmi_epsilon_one(xEnd,alpha,beta)

!   Trick to return x such that rho_reh=rho_end

    x = gmssmi_x_star(alpha,beta,wrad,junk,Pstar)    
    potStar = gmssmi_norm_potential(x,alpha,beta)
    epsOneStar = gmssmi_epsilon_one(x,alpha,beta)


    
    if (.not.slowroll_validity(epsOneStar)) stop 'gmssmi_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    gmssmi_lnrhoend = lnRhoEnd

  end function gmssmi_lnrhoend

  
end module gmssmireheat

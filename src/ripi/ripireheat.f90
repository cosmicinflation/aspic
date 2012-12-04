!renormalizable inflection point reheating functions in the slow-roll approximations

module ripireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use ripisr, only : ripi_epsilon_one, ripi_epsilon_two, ripi_epsilon_three
  use ripisr, only : ripi_norm_potential
  use ripisr, only : ripi_x_endinf, ripi_efold_primitive, ripi_x_trajectory
  implicit none

  private

  public ripi_x_star, ripi_lnrhoend 

contains

!returns x=phi/alpha such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function ripi_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: ripi_x_star
    real(kp), intent(in) :: alpha,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: ripiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = ripi_x_endinf(alpha)

    epsOneEnd = ripi_epsilon_one(xEnd,alpha)
    potEnd = ripi_norm_potential(xEnd,alpha)

    primEnd = ripi_efold_primitive(xEnd,alpha)

   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    ripiData%real1 = alpha
    ripiData%real2 = w
    ripiData%real3 = calF + primEnd

    mini = xEnd*(1._kp+epsilon(1._kp))
    maxi = 4._kp/(3._kp*alpha)*(1._kp-epsilon(1._kp)) !Position of the flat inflection point


    x = zbrent(find_ripi_x_star,mini,maxi,tolzbrent,ripiData)
    ripi_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (ripi_efold_primitive(x,alpha) - primEnd)
    endif

  end function ripi_x_star

  function find_ripi_x_star(x,ripiData)   
    implicit none
    real(kp) :: find_ripi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: ripiData

    real(kp) :: primStar,alpha,w,CalFplusprimEnd,potStar,epsOneStar

    alpha = ripiData%real1
    w = ripiData%real2
    CalFplusprimEnd = ripiData%real3

    primStar = ripi_efold_primitive(x,alpha)
    epsOneStar = ripi_epsilon_one(x,alpha)
    potStar = ripi_norm_potential(x,alpha)

    find_ripi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_ripi_x_star



  function ripi_lnrhoend(alpha,Pstar) 
    implicit none
    real(kp) :: ripi_lnrhoend
    real(kp), intent(in) :: alpha,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = ripi_x_endinf(alpha)

    potEnd  = ripi_norm_potential(xEnd,alpha)

    epsOneEnd = ripi_epsilon_one(xEnd,alpha)

!   Trick to return x such that rho_reh=rho_end

    x = ripi_x_star(alpha,wrad,junk,Pstar)  

    potStar = ripi_norm_potential(x,alpha)
    epsOneStar = ripi_epsilon_one(x,alpha)

    
    if (.not.slowroll_validity(epsOneStar)) stop 'ripi_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    ripi_lnrhoend = lnRhoEnd

  end function ripi_lnrhoend

  
end module ripireheat

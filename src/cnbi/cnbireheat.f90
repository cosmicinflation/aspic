!constant ns B reheating functions in the slow-roll approximations

module cnbireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use cnbisr, only : cnbi_epsilon_one, cnbi_epsilon_two, cnbi_epsilon_three
  use cnbisr, only : cnbi_norm_potential, cnbi_x_max
  use cnbisr, only : cnbi_x_endinf, cnbi_efold_primitive, cnbi_x_trajectory
  implicit none

  private

  public cnbi_x_star, cnbi_lnrhoend 

contains

!returns x=phi/alpha such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function cnbi_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: cnbi_x_star
    real(kp), intent(in) :: alpha,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: cnbiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = cnbi_x_endinf(alpha)

    epsOneEnd = cnbi_epsilon_one(xEnd,alpha)
    potEnd = cnbi_norm_potential(xEnd,alpha)

    primEnd = cnbi_efold_primitive(xEnd,alpha)

    
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)


    cnbiData%real1 = alpha
    cnbiData%real2 = w
    cnbiData%real3 = calF + primEnd

    mini = xEnd*(1._kp+epsilon(1._kp))
    maxi = cnbi_x_max(alpha)*(1._kp-epsilon(1._kp)) !Position where slow roll stops being valid


    x = zbrent(find_cnbi_x_star,mini,maxi,tolzbrent,cnbiData)
    cnbi_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (cnbi_efold_primitive(x,alpha) - primEnd)
    endif

  end function cnbi_x_star

  function find_cnbi_x_star(x,cnbiData)   
    implicit none
    real(kp) :: find_cnbi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: cnbiData

    real(kp) :: primStar,alpha,w,CalFplusprimEnd,potStar,epsOneStar

    alpha = cnbiData%real1
    w = cnbiData%real2
    CalFplusprimEnd = cnbiData%real3

    primStar = cnbi_efold_primitive(x,alpha)
    epsOneStar = cnbi_epsilon_one(x,alpha)
    potStar = cnbi_norm_potential(x,alpha)

    find_cnbi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_cnbi_x_star



  function cnbi_lnrhoend(alpha,Pstar) 
    implicit none
    real(kp) :: cnbi_lnrhoend
    real(kp), intent(in) :: alpha,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = cnbi_x_endinf(alpha)

    potEnd  = cnbi_norm_potential(xEnd,alpha)

    epsOneEnd = cnbi_epsilon_one(xEnd,alpha)



!   Trick to return x such that rho_reh=rho_end

    x = cnbi_x_star(alpha,wrad,junk,Pstar)  

    potStar = cnbi_norm_potential(x,alpha)
    epsOneStar = cnbi_epsilon_one(x,alpha)

    
    if (.not.slowroll_validity(epsOneStar)) stop 'cnbi_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    cnbi_lnrhoend = lnRhoEnd

  end function cnbi_lnrhoend

  
end module cnbireheat

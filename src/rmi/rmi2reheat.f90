!running mass 2 reheating functions in the slow-roll approximations

module rmi2reheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use rmicomreh, only : rmi_x_star, rmi_x_rrad, rmi_x_rreh
  use rmi2sr, only : rmi2_epsilon_one, rmi2_epsilon_two, rmi2_epsilon_three
  use rmi2sr, only : rmi2_norm_potential, rmi2_efold_primitive

  implicit none

  private

  public rmi2_x_star, rmi2_lnrhoreh_max
  public rmi2_x_rrad, rmi2_x_rreh

contains



!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function rmi2_x_star(c,phi0,xend,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: rmi2_x_star
    real(kp), intent(in) :: c,phi0,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar
    
    real(kp) :: mini,maxi

    mini = 1._kp*(1._kp+epsilon(1._kp))
    maxi = xend*(1._kp-epsilon(1._kp))
   
    rmi2_x_star = rmi_x_star(c,phi0,xend,w,lnRhoReh,Pstar,mini,maxi,bfoldstar)

  end function rmi2_x_star


!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function rmi2_x_rrad(c,phi0,xend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: rmi2_x_rrad
    real(kp), intent(in) :: c,phi0,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: mini,maxi

    mini = 1._kp*(1._kp+epsilon(1._kp))
    maxi = xend*(1._kp-epsilon(1._kp))

    rmi2_x_rrad = rmi_x_rrad(c,phi0,xend,lnRrad,Pstar,mini,maxi,bfoldstar)

  end function rmi2_x_rrad

 
!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function rmi2_x_rreh(c,phi0,xend,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: rmi2_x_rreh
    real(kp), intent(in) :: c,phi0,xend,lnRreh
    real(kp), intent(out), optional :: bfoldstar
    
    real(kp) :: mini,maxi

    mini = 1._kp*(1._kp+epsilon(1._kp))
    maxi = xend*(1._kp-epsilon(1._kp))
   
    rmi2_x_rreh = rmi_x_rreh(c,phi0,xend,lnRreh,mini,maxi,bfoldstar)

  end function rmi2_x_rreh



  function rmi2_lnrhoreh_max(c,phi0,xend,Pstar) 
    implicit none
    real(kp) :: rmi2_lnrhoreh_max
    real(kp), intent(in) :: c, phi0, xend, Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    

    potEnd  = rmi2_norm_potential(xEnd,c,phi0)

    epsOneEnd = rmi2_epsilon_one(xEnd,c,phi0)


!   Trick to return x such that rho_reh=rho_end

    x = rmi2_x_star(c,phi0,xend,wrad,junk,Pstar) 

    potStar = rmi2_norm_potential(x,c,phi0)
    epsOneStar = rmi2_epsilon_one(x,c,phi0)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'rmi2_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    rmi2_lnrhoreh_max = lnRhoEnd

  end function rmi2_lnrhoreh_max

  
end module rmi2reheat

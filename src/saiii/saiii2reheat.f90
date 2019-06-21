!string axion I inflation reheating functions for inflation occuring
!at increasing field values x > xVmax in the slow-roll approximations

module saii2reheat
  use infprec, only : kp
  use srreheat, only : slowroll_validity
  use srreheat, only : ln_rho_endinf, ln_rho_reheat
  use saiicomreh, only : saii_x_star, saii_x_rrad, saii_x_rreh
  use saii2sr, only : saii2_epsilon_one, saii2_epsilon_two, saii2_epsilon_three
  use saii2sr, only : saii2_norm_potential, saii2_numacc_xinimin

  implicit none

  private

  public saii2_x_star, saii2_lnrhoreh_max
  public saii2_x_rrad, saii2_x_rreh

contains


!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function saii2_x_star(alpha,mu,xend,w,lnRhoReh,Pstar,bfoldstar)
    implicit none
    real(kp) :: saii2_x_star
    real(kp), intent(in) :: alpha,mu,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar
    
    real(kp) :: mini,maxi
   
    mini = saii2_numacc_xinimin(alpha,mu)
    maxi = xEnd

    saii2_x_star = saii_x_star(alpha,mu,w,lnRhoReh,Pstar,xend,mini,maxi,bfoldstar)

  end function saii2_x_star



!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function saii2_x_rrad(alpha,mu,xend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: saii2_x_rrad
    real(kp), intent(in) :: alpha,mu,xend,lnRrad,Pstar    
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: mini,maxi

    mini = saii2_numacc_xinimin(alpha,mu)
    maxi = xEnd
    
    saii2_x_rrad = saii_x_rrad(alpha,mu,lnRrad,Pstar,xend,mini,maxi,bfoldstar)

  end function saii2_x_rrad


!returns x given potential parameters, scalar power, and lnR.
!If present, returns the corresponding bfoldstar
  function saii2_x_rreh(alpha,mu,xend,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: saii2_x_rreh
    real(kp), intent(in) :: alpha,mu,xend,lnRreh 
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: mini,maxi

    mini = saii2_numacc_xinimin(alpha,mu)
    maxi = xEnd
    
    saii2_x_rreh = saii_x_rreh(alpha,mu,lnRreh,xend,mini,maxi,bfoldstar)

  end function saii2_x_rreh



  function saii2_lnrhoreh_max(alpha,mu,xend,Pstar) 
    implicit none
    real(kp) :: saii2_lnrhoreh_max
    real(kp), intent(in) :: alpha,mu,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd
    
    potEnd  = saii2_norm_potential(xEnd,alpha,mu)

    epsOneEnd = saii2_epsilon_one(xEnd,alpha,mu)

!   Trick to return x such that rho_reh=rho_end

    x = saii2_x_star(alpha,mu,xend,wrad,junk,Pstar)  

    potStar = saii2_norm_potential(x,alpha,mu)
    epsOneStar = saii2_epsilon_one(x,alpha,mu)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'saii2_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    saii2_lnrhoreh_max = lnRhoEnd

  end function saii2_lnrhoreh_max

  
end module saii2reheat

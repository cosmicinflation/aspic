!string axion I inflation reheating functions for inflation occuring
!at decreasing field values x < xVmax in the slow-roll approximations

module saii1reheat
  use infprec, only : kp
  use srreheat, only : slowroll_validity
  use srreheat, only : ln_rho_endinf, ln_rho_reheat
  use saiicomreh, only : saii_x_star, saii_x_rrad, saii_x_rreh
  use saii1sr, only : saii1_epsilon_one, saii1_epsilon_two, saii1_epsilon_three
  use saii1sr, only : saii1_norm_potential, saii1_numacc_xinimax

  implicit none

  private

  public saii1_x_star, saii1_lnrhoreh_max
  public saii1_x_rrad, saii1_x_rreh

contains


!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function saii1_x_star(alpha,mu,xend,w,lnRhoReh,Pstar,bfoldstar)
    implicit none
    real(kp) :: saii1_x_star
    real(kp), intent(in) :: alpha,mu,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar
    
    real(kp) :: mini,maxi
   
    mini = xEnd
    maxi = saii1_numacc_xinimax(alpha,mu)

    saii1_x_star = saii_x_star(alpha,mu,w,lnRhoReh,Pstar,xend,mini,maxi,bfoldstar)

  end function saii1_x_star



!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function saii1_x_rrad(alpha,mu,xend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: saii1_x_rrad
    real(kp), intent(in) :: alpha,mu,xend,lnRrad,Pstar    
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: mini,maxi

    mini = xEnd
    maxi = saii1_numacc_xinimax(alpha,mu)
    
    saii1_x_rrad = saii_x_rrad(alpha,mu,lnRrad,Pstar,xend,mini,maxi,bfoldstar)

  end function saii1_x_rrad


!returns x given potential parameters, scalar power, and lnR.
!If present, returns the corresponding bfoldstar
  function saii1_x_rreh(alpha,mu,xend,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: saii1_x_rreh
    real(kp), intent(in) :: alpha,mu,xend,lnRreh 
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: mini,maxi

    mini = xEnd
    maxi = saii1_numacc_xinimax(alpha,mu)
    
    saii1_x_rreh = saii_x_rreh(alpha,mu,lnRreh,xend,mini,maxi,bfoldstar)

  end function saii1_x_rreh



  function saii1_lnrhoreh_max(alpha,mu,xend,Pstar) 
    implicit none
    real(kp) :: saii1_lnrhoreh_max
    real(kp), intent(in) :: alpha,mu,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd
    
    potEnd  = saii1_norm_potential(xEnd,alpha,mu)

    epsOneEnd = saii1_epsilon_one(xEnd,alpha,mu)

!   Trick to return x such that rho_reh=rho_end

    x = saii1_x_star(alpha,mu,xend,wrad,junk,Pstar)  

    potStar = saii1_norm_potential(x,alpha,mu)
    epsOneStar = saii1_epsilon_one(x,alpha,mu)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'saii1_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    saii1_lnrhoreh_max = lnRhoEnd

  end function saii1_lnrhoreh_max

  
end module saii1reheat

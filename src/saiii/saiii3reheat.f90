!string axion II inflation reheating functions for inflation occuring
!at decreasing field values x < xVmax in the slow-roll approximations

module saiii3reheat
  use infprec, only : kp
  use srreheat, only : slowroll_validity
  use srreheat, only : ln_rho_endinf, ln_rho_reheat
  use saiiicomreh, only : saiii_x_star, saiii_x_rrad, saiii_x_rreh, saiiiXBig
  use saiii3sr, only : saiii3_epsilon_one, saiii3_epsilon_two, saiii3_epsilon_three
  use saiii3sr, only : saiii3_norm_potential

  implicit none

  private

  public saiii3_x_star, saiii3_lnrhoreh_max
  public saiii3_x_rrad, saiii3_x_rreh

  
contains


!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function saiii3_x_star(alpha,beta,mu,xend,w,lnRhoReh,Pstar,bfoldstar)
    implicit none
    real(kp) :: saiii3_x_star
    real(kp), intent(in) :: alpha,beta,mu,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar
    
    real(kp) :: mini,maxi
   
    mini = xEnd
    maxi = saiiiXBig

    saiii3_x_star = saiii_x_star(alpha,beta,mu,w,lnRhoReh,Pstar,xend,mini,maxi,bfoldstar)

  end function saiii3_x_star



!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function saiii3_x_rrad(alpha,beta,mu,xend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: saiii3_x_rrad
    real(kp), intent(in) :: alpha,beta,mu,xend,lnRrad,Pstar    
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: mini,maxi

    mini = xEnd
    maxi = saiiiXBig
    
    saiii3_x_rrad = saiii_x_rrad(alpha,beta,mu,lnRrad,Pstar,xend,mini,maxi,bfoldstar)

  end function saiii3_x_rrad


!returns x given potential parameters, scalar power, and lnR.
!If present, returns the corresponding bfoldstar
  function saiii3_x_rreh(alpha,beta,mu,xend,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: saiii3_x_rreh
    real(kp), intent(in) :: alpha,beta,mu,xend,lnRreh 
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: mini,maxi

    mini = xEnd
    maxi = saiiiXBig
    
    saiii3_x_rreh = saiii_x_rreh(alpha,beta,mu,lnRreh,xend,mini,maxi,bfoldstar)

  end function saiii3_x_rreh



  function saiii3_lnrhoreh_max(alpha,beta,mu,xend,Pstar) 
    implicit none
    real(kp) :: saiii3_lnrhoreh_max
    real(kp), intent(in) :: alpha,beta,mu,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd
    
    potEnd  = saiii3_norm_potential(xEnd,alpha,beta,mu)

    epsOneEnd = saiii3_epsilon_one(xEnd,alpha,beta,mu)

!   Trick to return x such that rho_reh=rho_end

    x = saiii3_x_star(alpha,beta,mu,xend,wrad,junk,Pstar)  

    potStar = saiii3_norm_potential(x,alpha,beta,mu)
    epsOneStar = saiii3_epsilon_one(x,alpha,beta,mu)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'saiii3_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    saiii3_lnrhoreh_max = lnRhoEnd

  end function saiii3_lnrhoreh_max

  
end module saiii3reheat

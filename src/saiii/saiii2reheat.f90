!string axion II inflation reheating functions for inflation occuring
!at increasing field values x > xVmax and x < xVmin with V(xVmin) <=
!0, in the slow-roll approximations

module saiii2reheat
  use infprec, only : kp
  use srreheat, only : slowroll_validity
  use srreheat, only : ln_rho_endinf, ln_rho_reheat
  use saiiicomreh, only : saiii_x_star, saiii_x_rrad, saiii_x_rreh
  use saiii2sr, only : saiii2_epsilon_one, saiii2_epsilon_two, saiii2_epsilon_three
  use saiii2sr, only : saiii2_norm_potential, saiii2_numacc_xinimin

  implicit none

  private

  public saiii2_x_star, saiii2_lnrhoreh_max
  public saiii2_x_rrad, saiii2_x_rreh

contains


!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function saiii2_x_star(alpha,beta,mu,xend,w,lnRhoReh,Pstar,bfoldstar)
    implicit none
    real(kp) :: saiii2_x_star
    real(kp), intent(in) :: alpha,beta,mu,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar
    
    real(kp) :: mini,maxi
   
    mini = saiii2_numacc_xinimin(alpha,beta,mu)
    maxi = xEnd

    saiii2_x_star = saiii_x_star(alpha,beta,mu,w,lnRhoReh,Pstar,xend,mini,maxi,bfoldstar)

  end function saiii2_x_star



!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function saiii2_x_rrad(alpha,beta,mu,xend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: saiii2_x_rrad
    real(kp), intent(in) :: alpha,beta,mu,xend,lnRrad,Pstar    
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: mini,maxi

    mini = saiii2_numacc_xinimin(alpha,beta,mu)
    maxi = xEnd
    
    saiii2_x_rrad = saiii_x_rrad(alpha,beta,mu,lnRrad,Pstar,xend,mini,maxi,bfoldstar)

  end function saiii2_x_rrad


!returns x given potential parameters, scalar power, and lnR.
!If present, returns the corresponding bfoldstar
  function saiii2_x_rreh(alpha,beta,mu,xend,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: saiii2_x_rreh
    real(kp), intent(in) :: alpha,beta,mu,xend,lnRreh 
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: mini,maxi

    mini = saiii2_numacc_xinimin(alpha,beta,mu)
    maxi = xEnd
    
    saiii2_x_rreh = saiii_x_rreh(alpha,beta,mu,lnRreh,xend,mini,maxi,bfoldstar)

  end function saiii2_x_rreh



  function saiii2_lnrhoreh_max(alpha,beta,mu,xend,Pstar) 
    implicit none
    real(kp) :: saiii2_lnrhoreh_max
    real(kp), intent(in) :: alpha,beta,mu,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd
    
    potEnd  = saiii2_norm_potential(xEnd,alpha,beta,mu)

    epsOneEnd = saiii2_epsilon_one(xEnd,alpha,beta,mu)

!   Trick to return x such that rho_reh=rho_end

    x = saiii2_x_star(alpha,beta,mu,xend,wrad,junk,Pstar)  

    potStar = saiii2_norm_potential(x,alpha,beta,mu)
    epsOneStar = saiii2_epsilon_one(x,alpha,beta,mu)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'saiii2_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    saiii2_lnrhoreh_max = lnRhoEnd

  end function saiii2_lnrhoreh_max

  
end module saiii2reheat

!running mass 4 reheating functions in the slow-roll approximations

module rmi4reheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use rmicomreh, only : rmi_x_star, rmi_x_rrad, rmi_x_rreh
  use rmi4sr, only : rmi4_epsilon_one, rmi4_epsilon_two, rmi4_epsilon_three
  use rmi4sr, only : rmi4_norm_potential, rmi4_efold_primitive
  use rmi4sr, only : rmi4_x_epsonemax
  

  implicit none

  private


  public rmi4_x_star, rmi4_x_rrad, rmi4_x_rreh, rmi4_lnrhoreh_max

contains


!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function rmi4_x_star(c,phi0,xend,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: rmi4_x_star
    real(kp), intent(in) :: c,phi0,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar
    
    real(kp) :: mini,maxi

    mini = xend +epsilon(1._kp)
    maxi = rmi4_x_epsonemax(c,phi0)
    
    rmi4_x_star = rmi_x_star(c,phi0,xend,w,lnRhoReh,Pstar,mini,maxi,bfoldstar)

  end function rmi4_x_star


!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function rmi4_x_rrad(c,phi0,xend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: rmi4_x_rrad
    real(kp), intent(in) :: c,phi0,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: mini,maxi

    mini = xend +epsilon(1._kp)
    maxi =  rmi4_x_epsonemax(c,phi0)
   
    rmi4_x_rrad = rmi_x_rrad(c,phi0,xend,lnRrad,Pstar,mini,maxi,bfoldstar)

  end function rmi4_x_rrad

 
!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function rmi4_x_rreh(c,phi0,xend,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: rmi4_x_rreh
    real(kp), intent(in) :: c,phi0,xend,lnRreh
    real(kp), intent(out), optional :: bfoldstar
    
    real(kp) :: mini,maxi

    mini = xend + epsilon(1._kp)
    maxi =  rmi4_x_epsonemax(c,phi0)

   
    rmi4_x_rreh = rmi_x_rreh(c,phi0,xend,lnRreh,mini,maxi,bfoldstar)

  end function rmi4_x_rreh



  function rmi4_lnrhoreh_max(c,phi0,xend,Pstar) 
    implicit none
    real(kp) :: rmi4_lnrhoreh_max
    real(kp), intent(in) :: c, phi0, xend, Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd
    

    potEnd  = rmi4_norm_potential(xEnd,c,phi0)

    epsOneEnd = rmi4_epsilon_one(xEnd,c,phi0)


!   Trick to return x such that rho_reh=rho_end

    x = rmi4_x_star(c,phi0,xend,wrad,junk,Pstar)  


    potStar = rmi4_norm_potential(x,c,phi0)
    epsOneStar = rmi4_epsilon_one(x,c,phi0)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'rmi4_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    rmi4_lnrhoreh_max = lnRhoEnd

  end function rmi4_lnrhoreh_max

  
end module rmi4reheat

!R+R^2+alpha R^3 inflation reheating functions for alphamin < alpha < 0
!in the slow-roll approximations

module ccsi3reheat
  use infprec, only : kp
  use srreheat, only : slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use ccsicomreh, only : ccsi_x_star, ccsi_x_rrad, ccsi_x_rreh
  use ccsi3sr, only : ccsi3_epsilon_one, ccsi3_epsilon_two, ccsi3_epsilon_three
  use ccsi3sr, only : ccsi3_norm_potential
  use ccsi3sr, only : ccsi3_x_endinf, ccsi3_efold_primitive
  use ccsi3sr, only : ccsi3_xinimax, ccsi3_check_params

  implicit none

  private

  public ccsi3_x_star, ccsi3_lnrhoreh_max
  public ccsi3_x_rrad, ccsi3_x_rreh

contains


!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function ccsi3_x_star(alpha,xend,w,lnRhoReh,Pstar,bfoldstar)
    implicit none
    real(kp) :: ccsi3_x_star
    real(kp), intent(in) :: alpha,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar
    
    real(kp) :: mini,maxi

    if (.not.ccsi3_check_params(alpha)) then
       stop 'ccsi3_x_star: ccsi3 requires alpha>=0'
    endif
   
    mini = xEnd
    maxi = ccsi3_xinimax(alpha)


    ccsi3_x_star = ccsi_x_star(alpha,w,lnRhoReh,Pstar,xend,mini,maxi,bfoldstar)

  end function ccsi3_x_star



!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function ccsi3_x_rrad(alpha,xend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: ccsi3_x_rrad
    real(kp), intent(in) :: alpha,xend,lnRrad,Pstar    
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: mini,maxi

    if (.not.ccsi3_check_params(alpha)) then
       stop 'ccsi3_x_rrad: ccsi3 requires alpha>=0'
    endif

    mini = xEnd
    maxi = ccsi3_xinimax(alpha)
    
    ccsi3_x_rrad = ccsi_x_rrad(alpha,lnRrad,Pstar,xend,mini,maxi,bfoldstar)

  end function ccsi3_x_rrad


!returns x given potential parameters, scalar power, and lnR.
!If present, returns the corresponding bfoldstar
  function ccsi3_x_rreh(alpha,xend,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: ccsi3_x_rreh
    real(kp), intent(in) :: alpha,xend,lnRreh 
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: mini,maxi

    if (.not.ccsi3_check_params(alpha)) then
       stop 'ccsi3_x_rreh: ccsi3 requires alpha>=0'
    endif

    mini = xEnd
    maxi = ccsi3_xinimax(alpha)
    
    ccsi3_x_rreh = ccsi_x_rreh(alpha,lnRreh,xend,mini,maxi,bfoldstar)

  end function ccsi3_x_rreh



  function ccsi3_lnrhoreh_max(alpha,xend,Pstar) 
    implicit none
    real(kp) :: ccsi3_lnrhoreh_max
    real(kp), intent(in) :: alpha,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd
    
    potEnd  = ccsi3_norm_potential(xEnd,alpha)

    epsOneEnd = ccsi3_epsilon_one(xEnd,alpha)

!   Trick to return x such that rho_reh=rho_end

    x = ccsi3_x_star(alpha,xend,wrad,junk,Pstar)  

    potStar = ccsi3_norm_potential(x,alpha)
    epsOneStar = ccsi3_epsilon_one(x,alpha)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'ccsi3_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    ccsi3_lnrhoreh_max = lnRhoEnd

  end function ccsi3_lnrhoreh_max

  
end module ccsi3reheat

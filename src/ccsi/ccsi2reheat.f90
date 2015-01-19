!R+R^2+alpha R^3 inflation reheating functions for alpha>0 and x>xVmax
!in the slow-roll approximations

module ccsi2reheat
  use infprec, only : kp
  use srreheat, only : slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use ccsicomreh, only : ccsi_x_star, ccsi_x_rrad, ccsi_x_rreh
  use ccsi2sr, only : ccsi2_epsilon_one, ccsi2_epsilon_two, ccsi2_epsilon_three
  use ccsi2sr, only : ccsi2_norm_potential, ccsi2_efold_primitive
  use ccsi2sr, only : ccsi2_check_params, ccsi2_numacc_xinimin

  implicit none

  private

  public ccsi2_x_star, ccsi2_lnrhoreh_max
  public ccsi2_x_rrad, ccsi2_x_rreh

contains


!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function ccsi2_x_star(alpha,xend,w,lnRhoReh,Pstar,bfoldstar)
    implicit none
    real(kp) :: ccsi2_x_star
    real(kp), intent(in) :: alpha,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar
    
    real(kp) :: mini,maxi

    if (.not.ccsi2_check_params(alpha)) then
       stop 'ccsi2_x_star: ccsi2 requires alpha>=0'
    endif

    mini = ccsi2_numacc_xinimin(alpha)
    maxi = xend

    ccsi2_x_star = ccsi_x_star(alpha,w,lnRhoReh,Pstar,xend,mini,maxi,bfoldstar)

  end function ccsi2_x_star



!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function ccsi2_x_rrad(alpha,xend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: ccsi2_x_rrad
    real(kp), intent(in) :: alpha,xend,lnRrad,Pstar    
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: mini,maxi

    if (.not.ccsi2_check_params(alpha)) then
       stop 'ccsi2_x_rrad: ccsi2 requires alpha>=0'
    endif
  
    mini = ccsi2_numacc_xinimin(alpha)
    maxi = xend
   
    ccsi2_x_rrad = ccsi_x_rrad(alpha,lnRrad,Pstar,xend,mini,maxi,bfoldstar)

  end function ccsi2_x_rrad


!returns x given potential parameters, scalar power, and lnR.
!If present, returns the corresponding bfoldstar
  function ccsi2_x_rreh(alpha,xend,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: ccsi2_x_rreh
    real(kp), intent(in) :: alpha,xend,lnRreh 
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: mini,maxi

    if (.not.ccsi2_check_params(alpha)) then
       stop 'ccsi2_x_rreh: ccsi2 requires alpha>=0'
    endif
   
    mini = ccsi2_numacc_xinimin(alpha)
    maxi = xend
    
    ccsi2_x_rreh = ccsi_x_rreh(alpha,lnRreh,xend,mini,maxi,bfoldstar)

  end function ccsi2_x_rreh



  function ccsi2_lnrhoreh_max(alpha,xend,Pstar) 
    implicit none
    real(kp) :: ccsi2_lnrhoreh_max
    real(kp), intent(in) :: alpha,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd
    
    potEnd  = ccsi2_norm_potential(xEnd,alpha)

    epsOneEnd = ccsi2_epsilon_one(xEnd,alpha)

!   Trick to return x such that rho_reh=rho_end

    x = ccsi2_x_star(alpha,xend,wrad,junk,Pstar)  

    potStar = ccsi2_norm_potential(x,alpha)
    epsOneStar = ccsi2_epsilon_one(x,alpha)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'ccsi2_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    ccsi2_lnrhoreh_max = lnRhoEnd

  end function ccsi2_lnrhoreh_max

  
end module ccsi2reheat

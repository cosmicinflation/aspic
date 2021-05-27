!R+R^2+alpha R^3 inflation reheating functions for alpha>0 and x<xVmax
!in the slow-roll approximations

module ccsi1reheat
  use infprec, only : kp
  use srreheat, only : slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use ccsicomreh, only : ccsi_x_star, ccsi_x_rrad, ccsi_x_rreh
  use ccsi1sr, only : ccsi1_epsilon_one, ccsi1_epsilon_two, ccsi1_epsilon_three
  use ccsi1sr, only : ccsi1_norm_potential
  use ccsi1sr, only : ccsi1_x_endinf, ccsi1_efold_primitive
  use ccsi1sr, only : ccsi1_numacc_xinimax, ccsi1_check_params

  implicit none

  private

  public ccsi1_x_star, ccsi1_lnrhoreh_max
  public ccsi1_x_rrad, ccsi1_x_rreh

contains


!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function ccsi1_x_star(alpha,xend,w,lnRhoReh,Pstar,bfoldstar)
    implicit none
    real(kp) :: ccsi1_x_star
    real(kp), intent(in) :: alpha,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar
    
    real(kp) :: mini,maxi

    if (.not.ccsi1_check_params(alpha)) then
       stop 'ccsi1_x_star: ccsi1 requires alpha>=0'
    endif
   
    mini = xEnd
    maxi = ccsi1_numacc_xinimax(alpha)


    ccsi1_x_star = ccsi_x_star(alpha,w,lnRhoReh,Pstar,xend,mini,maxi,bfoldstar)

  end function ccsi1_x_star



!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function ccsi1_x_rrad(alpha,xend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: ccsi1_x_rrad
    real(kp), intent(in) :: alpha,xend,lnRrad,Pstar    
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: mini,maxi

    if (.not.ccsi1_check_params(alpha)) then
       stop 'ccsi1_x_rrad: ccsi1 requires alpha>=0'
    endif

    mini = xEnd
    maxi = ccsi1_numacc_xinimax(alpha)
    
    ccsi1_x_rrad = ccsi_x_rrad(alpha,lnRrad,Pstar,xend,mini,maxi,bfoldstar)

  end function ccsi1_x_rrad


!returns x given potential parameters, scalar power, and lnR.
!If present, returns the corresponding bfoldstar
  function ccsi1_x_rreh(alpha,xend,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: ccsi1_x_rreh
    real(kp), intent(in) :: alpha,xend,lnRreh 
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: mini,maxi

    if (.not.ccsi1_check_params(alpha)) then
       stop 'ccsi1_x_rreh: ccsi1 requires alpha>=0'
    endif

    mini = xEnd
    maxi = ccsi1_numacc_xinimax(alpha)
    
    ccsi1_x_rreh = ccsi_x_rreh(alpha,lnRreh,xend,mini,maxi,bfoldstar)

  end function ccsi1_x_rreh



  function ccsi1_lnrhoreh_max(alpha,xend,Pstar) 
    implicit none
    real(kp) :: ccsi1_lnrhoreh_max
    real(kp), intent(in) :: alpha,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd,lnOmega4End
    
    potEnd  = ccsi1_norm_potential(xEnd,alpha)

    epsOneEnd = ccsi1_epsilon_one(xEnd,alpha)
    
!   Trick to return x such that rho_reh=rho_end

    x = ccsi1_x_star(alpha,xend,wrad,junk,Pstar)  

    potStar = ccsi1_norm_potential(x,alpha)
    epsOneStar = ccsi1_epsilon_one(x,alpha)
    lnOmega4End = 2._kp*xend
    
    if (.not.slowroll_validity(epsOneStar)) stop 'ccsi1_lnrhoreh_max: slow-roll violated!'

!Jordan frame    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar,lnOmega4End)

    ccsi1_lnrhoreh_max = lnRhoEnd

  end function ccsi1_lnrhoreh_max

  
end module ccsi1reheat

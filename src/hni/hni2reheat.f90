!hybrid natural 2 reheating functions in the slow-roll approximation

module hni2reheat
  use infprec, only : kp
  use srreheat, only : slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf,ln_rho_reheat
  use hnicomreh, only : hni_x_star, hni_x_rrad, hni_x_rreh
  use hni2sr, only : hni2_epsilon_one, hni2_epsilon_two
  use hni2sr, only : hni2_norm_potential, hni2_check_params
  
  implicit none

  private

  public hni2_x_star, hni2_lnrhoreh_max
  public hni2_x_rrad, hni2_x_rreh

contains

!returns x given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the correspoding bfoldstar
  function hni2_x_star(alpha,f,xend,w,lnRhoReh,Pstar,bfold)    
    implicit none
    real(kp) :: hni2_x_star
    real(kp), intent(in) :: alpha,f,xend,w,lnRhoReh,Pstar
    real(kp), intent(out), optional :: bfold
   
    real(kp) :: mini,maxi
    
    if (.not.hni2_check_params(alpha,f)) then
       stop 'hni2_x_star: hni2 requires alpha < alphamin!'
    endif
   
    mini = epsilon(1._kp)
    maxi = xend

    hni2_x_star = hni_x_star(alpha,f,w,lnRhoReh,Pstar,xend,mini,maxi,bfold)  
    
    
  end function hni2_x_star


!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function hni2_x_rrad(alpha,f,xend,lnRrad,Pstar,bfold)    
    implicit none
    real(kp) :: hni2_x_rrad
    real(kp), intent(in) :: alpha,f,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfold
   
    real(kp) :: mini,maxi
    
    if (.not.hni2_check_params(alpha,f)) then
       stop 'hni2_x_rrad: hni2 requires alpha < alphamin!'
    endif

    mini = epsilon(1._kp)
    maxi = xend
    
    hni2_x_rrad = hni_x_rrad(alpha,f,lnRrad,Pstar,xend,mini,maxi,bfold) 
    
  end function hni2_x_rrad

  


!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function hni2_x_rreh(alpha,f,xend,lnRreh,bfold)    
    implicit none
    real(kp) :: hni2_x_rreh
    real(kp), intent(in) :: alpha,f,xend,lnRreh
    real(kp), intent(out), optional :: bfold

    real(kp) :: mini,maxi

     
    mini = epsilon(1._kp)
    maxi = xend
    
    hni2_x_rreh =  hni_x_rreh(alpha,f,lnRreh,xend,mini,maxi,bfold)

  end function hni2_x_rreh
 


  function hni2_lnrhoreh_max(alpha,f,xend,Pstar) 
    implicit none
    real(kp) :: hni2_lnrhoreh_max
    real(kp), intent(in) :: alpha,f,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp), parameter :: wrad = 1._kp/3._kp
    real(kp), parameter :: junk= 0._kp
    real(kp) :: lnRhoEnd
       
    potEnd  = hni2_norm_potential(xEnd,alpha,f)
    epsOneEnd = hni2_epsilon_one(xEnd,alpha,f)
       
    x = hni2_x_star(alpha,f,xend,wrad,junk,Pstar)

    potStar = hni2_norm_potential(x,alpha,f)
    epsOneStar = hni2_epsilon_one(x,alpha,f)
    
    if (.not.slowroll_validity(epsOneStar)) then
        print*,'xstar=',x,'  epsOneStar=',epsOneStar 
        stop 'hni2_lnrhoreh_max: slow-roll violated!'
    endif
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    hni2_lnrhoreh_max = lnRhoEnd

  end function hni2_lnrhoreh_max

  
  
end module hni2reheat

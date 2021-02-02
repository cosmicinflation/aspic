!hybrid natural 1 reheating functions in the slow-roll approximation

module hni1reheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf,ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use hnicomreh, only : hni_x_star, hni_x_rrad, hni_x_rreh
  use hni1sr, only : hni1_epsilon_one, hni1_epsilon_two, hni1_epsilon_three
  use hni1sr, only : hni1_norm_potential,hni1_efold_primitive,hni1_x_endinf
  implicit none

  private

  public hni1_x_star, hni1_lnrhoreh_max
  public hni1_x_rrad, hni1_x_rreh

contains

!returns x given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the correspoding bfoldstar
  function hni1_x_star(alpha,f,xend,w,lnRhoReh,Pstar,bfold)    
    implicit none
    real(kp) :: hni1_x_star
    real(kp), intent(in) :: alpha,f,xend,w,lnRhoReh,Pstar
    real(kp), intent(out), optional :: bfold
   
    real(kp) :: mini,maxi
    
    if (.not.hni1_check_params(alpha,f)) then
       stop 'hni1_x_star: hni1 requires alpha > alphamin!'
    endif

    mini = epsilon(1._kp)
    maxi = xend*(1._kp-epsilon(1._kp))

    hni1_x_star = hni_x_star(alpha,f,w,lnRhoReh,Pstar,xend,mini,maxi,bfold)  
    
    
  end function hni1_x_star


!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function hni1_x_rrad(alpha,f,xend,lnRrad,Pstar,bfold)    
    implicit none
    real(kp) :: hni1_x_rrad
    real(kp), intent(in) :: alpha,f,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfold
   
    real(kp) :: mini,maxi
    
    if (.not.hni1_check_params(alpha,f)) then
       stop 'hni1_x_rrad: hni1 requires alpha > alphamin!'
    endif

    mini = epsilon(1._kp)
    maxi = xend*(1._kp-epsilon(1._kp))
    
    hni1_x_rrad = hni_x_rrad(alpha,f,lnRrad,Pstar,xend,mini,maxi,bfold) 
    
  end function hni1_x_rrad

  


!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function hni1_x_rreh(alpha,f,xend,lnRreh,bfold)    
    implicit none
    real(kp) :: hni1_x_rreh
    real(kp), intent(in) :: alpha,f,xend,lnRreh
    real(kp), intent(out), optional :: bfold

    real(kp) :: mini,maxi

     if (.not.hni1_check_params(alpha,f)) then
       stop 'hni1_x_rreh: hni1 requires alpha > alphamin!'
    endif

    mini = epsilon(1._kp)
    maxi = xend*(1._kp-epsilon(1._kp))
    
    hni1_x_rreh =  hni_x_rreh(alpha,f,lnRreh,xend,mini,maxi,bfold)

  end function hni1_x_rreh
 


  function hni1_lnrhoreh_max(alpha,f,xend,Pstar) 
    implicit none
    real(kp) :: hni1_lnrhoreh_max
    real(kp), intent(in) :: alpha,f,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp), parameter :: wrad = 1._kp/3._kp
    real(kp), parameter :: junk= 0._kp
    real(kp) :: lnRhoEnd
        
    potEnd  = hni1_norm_potential(xEnd,alpha,f)
    epsOneEnd = hni1_epsilon_one(xEnd,alpha,f)
       
    x = hni1_x_star(alpha,f,xend,wrad,junk,Pstar)

    potStar = hni1_norm_potential(x,alpha,f)
    epsOneStar = hni1_epsilon_one(x,alpha,f)
    
    if (.not.slowroll_validity(epsOneStar)) then
        print*,'xstar=',x,'  epsOneStar=',epsOneStar 
        stop 'hni1_lnrhoreh_max: slow-roll violated!'
    endif
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    hni1_lnrhoreh_max = lnRhoEnd

  end function hni1_lnrhoreh_max

  
  
end module hni1reheat

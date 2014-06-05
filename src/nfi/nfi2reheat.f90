!N-formalism inflation 2 reheating functions in the slow-roll approximations

module nfi2reheat
  use infprec, only : kp
  use srreheat, only : slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use nficomreh, only : nfi_x_star, nfi_x_rrad, nfi_x_rreh
  use nfi2sr, only : nfi2_epsilon_one, nfi2_epsilon_two, nfi2_epsilon_three
  use nfi2sr, only : nfi2_norm_potential, nfi2_efold_primitive
  use nfi2sr, only : nfi2_check_params, nfi2_numacc_xinimax

  implicit none

  private

  public nfi2_x_star, nfi2_lnrhoreh_max
  public nfi2_x_rrad, nfi2_x_rreh

contains


!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function nfi2_x_star(a,b,xend,w,lnRhoReh,Pstar,bfoldstar)
    implicit none
    real(kp) :: nfi2_x_star
    real(kp), intent(in) :: a,b,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar
    
    real(kp) :: mini,maxi

    if (.not.nfi2_check_params(a,b)) then
       stop 'nfi2_x_star: nfi2 requires a<0, b>1'
    endif
       
    mini = xend
    maxi = nfi2_numacc_xinimax(a,b)

    nfi2_x_star = nfi_x_star(a,b,w,lnRhoReh,Pstar,xend,mini,maxi,bfoldstar)

  end function nfi2_x_star



!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function nfi2_x_rrad(a,b,xend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: nfi2_x_rrad
    real(kp), intent(in) :: a,b,xend,lnRrad,Pstar    
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: mini,maxi

    if (.not.nfi2_check_params(a,b)) then
       stop 'nfi2_x_rrad: nfi2 requires a<0, b>1'
    endif
  
    mini = xend
    maxi = nfi2_numacc_xinimax(a,b)
   
    nfi2_x_rrad = nfi_x_rrad(a,b,lnRrad,Pstar,xend,mini,maxi,bfoldstar)

  end function nfi2_x_rrad


!returns x given potential parameters, scalar power, and lnR.
!If present, returns the corresponding bfoldstar
  function nfi2_x_rreh(a,b,xend,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: nfi2_x_rreh
    real(kp), intent(in) :: a,b,xend,lnRreh 
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: mini,maxi

    if (.not.nfi2_check_params(a,b)) then
       stop 'nfi2_x_rreh: nfi2 requires a<0, b>1'
    endif
   
    mini = xend
    maxi = nfi2_numacc_xinimax(a,b)
    
    nfi2_x_rreh = nfi_x_rreh(a,b,lnRreh,xend,mini,maxi,bfoldstar)

  end function nfi2_x_rreh



  function nfi2_lnrhoreh_max(a,b,xend,Pstar) 
    implicit none
    real(kp) :: nfi2_lnrhoreh_max
    real(kp), intent(in) :: a,b,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd
    
    potEnd  = nfi2_norm_potential(xEnd,a,b)

    epsOneEnd = nfi2_epsilon_one(xEnd,a,b)

!   Trick to return x such that rho_reh=rho_end

    x = nfi2_x_star(a,b,xend,wrad,junk,Pstar)  

    potStar = nfi2_norm_potential(x,a,b)
    epsOneStar = nfi2_epsilon_one(x,a,b)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'nfi2_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    nfi2_lnrhoreh_max = lnRhoEnd

  end function nfi2_lnrhoreh_max

  
end module nfi2reheat

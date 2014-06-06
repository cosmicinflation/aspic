!N-formalism inflation 1 reheating functions in the slow-roll approximations

module nfi3reheat
  use infprec, only : kp
  use srreheat, only : slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use nficomreh, only : nfi_x_star, nfi_x_rrad, nfi_x_rreh
  use nfi3sr, only : nfi3_epsilon_one, nfi3_epsilon_two, nfi3_epsilon_three
  use nfi3sr, only : nfi3_norm_potential
  use nfi3sr, only : nfi3_x_endinf, nfi3_efold_primitive
  use nfi3sr, only : nfi3_check_params, nfi3_numacc_xinimax

  implicit none

  private

  public nfi3_x_star, nfi3_lnrhoreh_max
  public nfi3_x_rrad, nfi3_x_rreh

contains


!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function nfi3_x_star(a,b,w,lnRhoReh,Pstar,bfoldstar)
    implicit none
    real(kp) :: nfi3_x_star
    real(kp), intent(in) :: a,b,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar
    
    real(kp) :: mini,maxi
    real(kp) :: xend

    if (.not.nfi3_check_params(a,b)) then
       stop 'nfi3_x_star: nfi3 requires a>0, b>1'
    endif
   
    xEnd = nfi3_x_endinf(a,b)
    mini = xEnd
    maxi = nfi3_numacc_xinimax(a,b)

    nfi3_x_star = nfi_x_star(a,b,w,lnRhoReh,Pstar,xend,mini,maxi,bfoldstar)

  end function nfi3_x_star



!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function nfi3_x_rrad(a,b,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: nfi3_x_rrad
    real(kp), intent(in) :: a,b,lnRrad,Pstar    
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: mini,maxi
    real(kp) :: xend

    if (.not.nfi3_check_params(a,b)) then
       stop 'nfi3_x_rrad: nfi3 requires a>0, b>1'
    endif

    xEnd = nfi3_x_endinf(a,b)
    mini = xEnd
    maxi = nfi3_numacc_xinimax(a,b)
    
    
    nfi3_x_rrad = nfi_x_rrad(a,b,lnRrad,Pstar,xend,mini,maxi,bfoldstar)

  end function nfi3_x_rrad


!returns x given potential parameters, scalar power, and lnR.
!If present, returns the corresponding bfoldstar
  function nfi3_x_rreh(a,b,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: nfi3_x_rreh
    real(kp), intent(in) :: a,b,lnRreh 
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: mini,maxi
    real(kp) :: xend

    if (.not.nfi3_check_params(a,b)) then
       stop 'nfi3_x_rreh: nfi3 requires a>0, b>1'
    endif

    xEnd = nfi3_x_endinf(a,b)
    mini = xEnd
    maxi = nfi3_numacc_xinimax(a,b)
    
    nfi3_x_rreh = nfi_x_rreh(a,b,lnRreh,xend,mini,maxi,bfoldstar)

  end function nfi3_x_rreh



  function nfi3_lnrhoreh_max(a,b,Pstar) 
    implicit none
    real(kp) :: nfi3_lnrhoreh_max
    real(kp), intent(in) :: a,b,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd
    
    xEnd = nfi3_x_endinf(a,b)

    potEnd  = nfi3_norm_potential(xEnd,a,b)

    epsOneEnd = nfi3_epsilon_one(xEnd,a,b)

!   Trick to return x such that rho_reh=rho_end

    x = nfi3_x_star(a,b,wrad,junk,Pstar)  

    potStar = nfi3_norm_potential(x,a,b)
    epsOneStar = nfi3_epsilon_one(x,a,b)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'nfi3_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    nfi3_lnrhoreh_max = lnRhoEnd

  end function nfi3_lnrhoreh_max

  
end module nfi3reheat

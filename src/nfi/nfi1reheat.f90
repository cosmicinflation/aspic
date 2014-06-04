!N-formalism inflation 1 reheating functions in the slow-roll approximations

module nfi1reheat
  use infprec, only : kp
  use srreheat, only : slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use nficomreh, only : nfi_x_star, nfi_x_rrad, nfi_x_rreh
  use nfi1sr, only : nfi1_epsilon_one, nfi1_epsilon_two, nfi1_epsilon_three
  use nfi1sr, only : nfi1_norm_potential
  use nfi1sr, only : nfi1_x_endinf, nfi1_efold_primitive
  use nfi1sr, only : nfi1_numacc_xinimin, nfi1_check_params

  implicit none

  private

  public nfi1_x_star, nfi1_lnrhoreh_max
  public nfi1_x_rrad, nfi1_x_rreh

contains


!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function nfi1_x_star(a,b,w,lnRhoReh,Pstar,bfoldstar)
    implicit none
    real(kp) :: nfi1_x_star
    real(kp), intent(in) :: a,b,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar
    
    real(kp) :: mini,maxi
    real(kp) :: xend

    if (.not.nfi1_check_params(a,b)) then
       stop 'nfi1_x_star: nfi1 requires a>0, b>1'
    endif
   
    xEnd = nfi1_x_endinf(a,b)
    mini = nfi1_numacc_xinimin(a,b)
    maxi = xEnd

    nfi1_x_star = nfi_x_star(a,b,w,lnRhoReh,Pstar,xend,mini,maxi,bfoldstar)

  end function nfi1_x_star



!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function nfi1_x_rrad(a,b,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: nfi1_x_rrad
    real(kp), intent(in) :: a,b,lnRrad,Pstar    
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: mini,maxi
    real(kp) :: xend

    if (.not.nfi1_check_params(a,b)) then
       stop 'nfi1_x_rrad: nfi1 requires a>0, b>1'
    endif

    xEnd = nfi1_x_endinf(a,b)
    mini = nfi1_numacc_xinimin(a,b)
    maxi = xEnd
    
    nfi1_x_rrad = nfi_x_rrad(a,b,lnRrad,Pstar,xend,mini,maxi,bfoldstar)

  end function nfi1_x_rrad


!returns x given potential parameters, scalar power, and lnR.
!If present, returns the corresponding bfoldstar
  function nfi1_x_rreh(a,b,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: nfi1_x_rreh
    real(kp), intent(in) :: a,b,lnRreh 
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: mini,maxi
    real(kp) :: xend

    if (.not.nfi1_check_params(a,b)) then
       stop 'nfi1_x_rreh: nfi1 requires a>0, b>1'
    endif

    xEnd = nfi1_x_endinf(a,b)
    mini = nfi1_numacc_xinimin(a,b)
    maxi = xEnd
    
    nfi1_x_rreh = nfi_x_rreh(a,b,lnRreh,xend,mini,maxi,bfoldstar)

  end function nfi1_x_rreh



  function nfi1_lnrhoreh_max(a,b,Pstar) 
    implicit none
    real(kp) :: nfi1_lnrhoreh_max
    real(kp), intent(in) :: a,b,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd
    
    xEnd = nfi1_x_endinf(a,b)

    potEnd  = nfi1_norm_potential(xEnd,a,b)

    epsOneEnd = nfi1_epsilon_one(xEnd,a,b)

!   Trick to return x such that rho_reh=rho_end

    x = nfi1_x_star(a,b,wrad,junk,Pstar)  

    potStar = nfi1_norm_potential(x,a,b)
    epsOneStar = nfi1_epsilon_one(x,a,b)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'nfi1_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    nfi1_lnrhoreh_max = lnRhoEnd

  end function nfi1_lnrhoreh_max

  
end module nfi1reheat

!N-formalism inflation 4 reheating functions in the slow-roll approximations

module nfi4reheat
  use infprec, only : kp
  use srreheat, only : slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use nficomreh, only : nfi_x_star, nfi_x_rrad, nfi_x_rreh
  use nfi4sr, only : nfi4_epsilon_one, nfi4_epsilon_two, nfi4_epsilon_three
  use nfi4sr, only : nfi4_norm_potential, nfi4_efold_primitive
  use nfi4sr, only : nfi4_check_params, nfi4_xinimin

  implicit none

  private

  public nfi4_x_star, nfi4_lnrhoreh_max
  public nfi4_x_rrad, nfi4_x_rreh

contains


!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function nfi4_x_star(a,b,xend,w,lnRhoReh,Pstar,bfoldstar)
    implicit none
    real(kp) :: nfi4_x_star
    real(kp), intent(in) :: a,b,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar
    
    real(kp) :: mini,maxi

    if (.not.nfi4_check_params(a,b)) then
       stop 'nfi4_x_star: nfi4 requires a>0, 0<b<1 or a<0, b<0'
    endif
       
    mini = nfi4_xinimin(a,b)
    maxi = xend

    nfi4_x_star = nfi_x_star(a,b,w,lnRhoReh,Pstar,xend,mini,maxi,bfoldstar)

  end function nfi4_x_star



!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function nfi4_x_rrad(a,b,xend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: nfi4_x_rrad
    real(kp), intent(in) :: a,b,xend,lnRrad,Pstar    
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: mini,maxi

    if (.not.nfi4_check_params(a,b)) then
       stop 'nfi4_x_rrad: nfi4 requires a>0, 0<b<1 or a<0, b<0'
    endif
  
    mini = nfi4_xinimin(a,b)
    maxi = xend
   
    nfi4_x_rrad = nfi_x_rrad(a,b,lnRrad,Pstar,xend,mini,maxi,bfoldstar)

  end function nfi4_x_rrad


!returns x given potential parameters, scalar power, and lnR.
!If present, returns the corresponding bfoldstar
  function nfi4_x_rreh(a,b,xend,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: nfi4_x_rreh
    real(kp), intent(in) :: a,b,xend,lnRreh 
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: mini,maxi

    if (.not.nfi4_check_params(a,b)) then
       stop 'nfi4_x_rreh: nfi4 requires a>0, 0<b<1 or a<0, b<0'
    endif
   
    mini = nfi4_xinimin(a,b)
    maxi = xend
    
    nfi4_x_rreh = nfi_x_rreh(a,b,lnRreh,xend,mini,maxi,bfoldstar)

  end function nfi4_x_rreh



  function nfi4_lnrhoreh_max(a,b,xend,Pstar) 
    implicit none
    real(kp) :: nfi4_lnrhoreh_max
    real(kp), intent(in) :: a,b,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd
    
    potEnd  = nfi4_norm_potential(xEnd,a,b)

    epsOneEnd = nfi4_epsilon_one(xEnd,a,b)

!   Trick to return x such that rho_reh=rho_end

    x = nfi4_x_star(a,b,xend,wrad,junk,Pstar)  

    potStar = nfi4_norm_potential(x,a,b)
    epsOneStar = nfi4_epsilon_one(x,a,b)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'nfi4_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    nfi4_lnrhoreh_max = lnRhoEnd

  end function nfi4_lnrhoreh_max

  
end module nfi4reheat

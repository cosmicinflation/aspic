!Reheating functions for Non-Minimal Large Field Inflation 2,
!gracefully ending in the large field region x > xVmax

module nmlfi2reheat
  use infprec, only :  kp
  use nmlficommon, only : nmlfi_x, nmlfi_hbar, hbarBig, hbarSmall
  use nmlficommon, only : nmlfi_hbar_potmax

  use nmlficomreh, only : nmlfi_hbar_rreh, nmlfi_hbar_star, nmlfi_hbar_rrad
  use nmlficomreh, only : nmlfi_lnrhoreh_max
  
  use nmlfi2sr, only : nmlfi2_check_params
  
  implicit none

  private

!take as input the parametric field hbar  
  public nmlfi2_parametric_lnrhoreh_max
  public nmlfi2_hbar_star, nmlfi2_hbar_rrad, nmlfi2_hbar_rreh

!for convenience only
  public nmlfi2_x_star, nmlfi2_x_rrad, nmlfi2_x_rreh, nmlfi2_lnrhoreh_max

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!the parametric reheating functions, the ones you should use !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function nmlfi2_hbar_star(xi,p,hbarend,w,lnRhoReh,Pstar,bfoldstar)
    implicit none
    real(kp) :: nmlfi2_hbar_star
    real(kp), intent(in) :: xi,p,hbarend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: hbarmin, hbarmax

    if (.not.nmlfi2_check_params(p,xi)) then
       write(*,*)'xi= p= ',xi,p
       stop 'nmlfi2_hbar_star: NMLFI2 does not exist for these parameters!'
    endif
    
    hbarmin = nmlfi_hbar_potmax(xi,p)
    hbarmax = hbarend
    
    nmlfi2_hbar_star = nmlfi_hbar_star(xi,p,w,lnRhoReh,Pstar,hbarend,hbarmin,hbarmax,bfoldstar)
    
  end function nmlfi2_hbar_star



  function nmlfi2_hbar_rrad(xi,p,hbarend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: nmlfi2_hbar_rrad
    real(kp), intent(in) :: xi,p,hbarend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: hbarmin, hbarmax

    if (.not.nmlfi2_check_params(p,xi)) then
       write(*,*)'xi= p= ',xi,p
       stop 'nmlfi2_hbar_rrad: NMLFI2 does not exist for these parameters!'
    endif

    hbarmin = nmlfi_hbar_potmax(xi,p)
    hbarmax = hbarend
    
    nmlfi2_hbar_rrad = nmlfi_hbar_rrad(xi,p,lnRrad,Pstar,hbarend,hbarmin,hbarmax,bfoldstar)

  end function nmlfi2_hbar_rrad
  


  function nmlfi2_hbar_rreh(xi,p,hbarend,lnRreh,bfoldstar)
    implicit none
    real(kp) :: nmlfi2_hbar_rreh
    real(kp), intent(in) :: xi,p,hbarend,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: hbarmin, hbarmax
    
    if (.not.nmlfi2_check_params(p,xi)) then
       write(*,*)'xi= p= ',xi,p
       stop 'nmlfi2_hbar_rreh: NMLFI2 does not exist for these parameters!'
    endif

    hbarmin = nmlfi_hbar_potmax(xi,p)
    hbarmax = hbarend
   
    nmlfi2_hbar_rreh = nmlfi_hbar_rreh(xi,p,lnRreh,hbarend,hbarmin,hbarmax,bfoldstar)
        
  end function nmlfi2_hbar_rreh

  
    

!jordan frame
  function nmlfi2_parametric_lnrhoreh_max(xi,p,hbarend,Pstar)
    implicit none
    real(kp) :: nmlfi2_parametric_lnrhoreh_max
    real(kp), intent(in) :: xi,p,hbarend,Pstar
    real(kp) :: hbarmin, hbarmax
        
    hbarmin = nmlfi_hbar_potmax(xi,p)
    hbarmax = hbarend
    
    nmlfi2_parametric_lnrhoreh_max = nmlfi_lnrhoreh_max(xi,p,Pstar,hbarend,hbarmin,hbarmax)

  end function nmlfi2_parametric_lnrhoreh_max


!Convenience functions: inefficient wrapper taking as input the
!canonical field. For computations it is better to call the parametric
!functions directly, they work in terms of hbar rather than x to
!prevent unnecessary numerical inversions
  function nmlfi2_x_star(xi,p,xend,w,lnRhoReh,Pstar,bfoldstar)
    implicit none
    real(kp) :: nmlfi2_x_star
    real(kp), intent(in) :: xi,p,xend,w,lnRhoReh,Pstar
    real(kp), intent(out), optional :: bfoldstar
    
    real(kp) :: hbarstar, hbarend

    hbarend = nmlfi_hbar(xend,xi)
    hbarstar = nmlfi2_hbar_star(xi,p,hbarend,w,lnRhoReh,Pstar,bfoldstar)

    nmlfi2_x_star = nmlfi_x(hbarstar,xi)
    
  end function nmlfi2_x_star


  function nmlfi2_x_rrad(xi,p,xend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: nmlfi2_x_rrad
    real(kp), intent(in) :: xi,p,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar
    
    real(kp) :: hbarstar,hbarend

    hbarend = nmlfi_hbar(xend,xi)
    hbarstar = nmlfi2_hbar_rrad(xi,p,hbarend,lnRrad,Pstar,bfoldstar)

    nmlfi2_x_rrad = nmlfi_x(hbarstar,xi)

  end function nmlfi2_x_rrad

    
  function nmlfi2_x_rreh(xi,p,xend,lnRreh,bfoldstar)
    implicit none
    real(kp) :: nmlfi2_x_rreh
    real(kp), intent(in) :: xi,p,xend,lnRreh
    real(kp), intent(out), optional :: bfoldstar
    
    real(kp) :: hbarstar,hbarend

    hbarend = nmlfi_hbar(xend,xi)
    hbarstar = nmlfi2_hbar_rreh(xi,p,hbarend,lnRreh,bfoldstar)

    nmlfi2_x_rreh = nmlfi_x(hbarstar,xi)
    
  end function nmlfi2_x_rreh
  

  function nmlfi2_lnrhoreh_max(xi,p,xend,Pstar)
    implicit none
    real(kp) :: nmlfi2_lnrhoreh_max
    real(kp), intent(in) :: xi,p,xend,Pstar
    real(kp) :: hbarmin, hbarmax, hbarend
    
    hbarend = nmlfi_hbar(xend,xi)
  
    nmlfi2_lnrhoreh_max = nmlfi2_parametric_lnrhoreh_max(xi,p,hbarend,Pstar)

  end function nmlfi2_lnrhoreh_max

  
end module nmlfi2reheat

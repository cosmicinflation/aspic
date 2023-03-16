!Reheating functions for Non-Minimal Large Field Inflation 3, in the
!large field region x > xVmax where xend is an extra model parameter

module nmlficomreh
  use infprec, only :  kp
  use nmlficommon, only : nmlfi_x, nmlfi_hbar, hbarBig, hbarSmall
  use nmlfi3sr, only : nmlfi3_check_params
  
  implicit none

  private

!take as input the parametric field hbar  
  public nmlfi3_lnrhoreh_max
  public nmlfi3_hbar_star, nmlfi3_hbar_rrad, nmlfi3_hbar_rreh

!for convenience only
  public nmlfi3_x_star, nmlfi3_x_rrad, nmlfi3_x_rreh

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!the parametric reheating functions, the ones you should use !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function nmlfi3_hbar_star(xi,p,hbarend,w,lnRhoReh,Pstar,bfoldstar)
    implicit none
    real(kp) :: nmlfi3_hbar_star
    real(kp), intent(in) :: xi,p,hbarend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: hbarmin, hbarmax

    if (.not.nmlfi3_check_params(p,xi)) then
       write(*,*)'xi= p= ',xi,p
       stop 'nmlfi3_hbar_star: NMLFI3 does not exist for these parameters!'
    endif
    
    hbarmin = nmlfi_hbar_potmax(xi,p)
    hbarmax = hbarend
    
    nmlfi3_hbar_star = nmlfi_hbar_star(xi,p,w,lnRhoReh,Pstar,hbarend,hbarmin,hbarmax,bfoldstar)
    
  end function nmlfi3_hbar_star



  function nmlfi3_hbar_rrad(xi,p,hbarend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: nmlfi3_hbar_rrad
    real(kp), intent(in) :: xi,p,hbarend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: hbarmin, hbarmax

    if (.not.nmlfi3_check_params(p,xi)) then
       write(*,*)'xi= p= ',xi,p
       stop 'nmlfi3_hbar_rrad: NMLFI3 does not exist for these parameters!'
    endif

    hbarmin = nmlfi_hbar_potmax(xi,p)
    hbarmax = hbarend
    
    nmlfi3_hbar_rrad = nmlfi_hbar_rrad(xi,p,lnRrad,Pstar,hbarend,hbarmin,hbarmax,bfoldstar)

  end function nmlfi3_hbar_rrad
  


  function nmlfi3_hbar_rreh(xi,p,hbarend,lnRreh,bfoldstar)
    implicit none
    real(kp) :: nmlfi3_hbar_rreh
    real(kp), intent(in) :: xi,p,hbarend,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: hbarmin, hbarmax
    
    if (.not.nmlfi3_check_params(p,xi)) then
       write(*,*)'xi= p= ',xi,p
       stop 'nmlfi3_hbar_rreh: NMLFI3 does not exist for these parameters!'
    endif

    hbarmin = nmlfi_hbar_potmax(xi,p)
    hbarmax = hbarend
   
    nmlfi3_hbar_rreh = nmlfi_hbar_rreh(xi,p,lnRreh,hbarend,hbarmin,hbarmax,bfoldstar)
        
  end function nmlfi3_hbar_rreh

  
    
!jordan frame
  function nmlfi3_lnrhoreh_max(xi,p,hbarend,Pstar)
    implicit none
    real(kp) :: nmlfi3_lnrhoreh_max
    real(kp), intent(in) :: xi,p,hbarend,Pstar

    nmlfi3_lnrhoreh_max = nmlfi_lnrhoreh_max(xi,p,hbarend,Pstar)
    
  end function nmlfi3_lnrhoreh_max



!Convenience functions: inefficient wrapper taking as input the
!canonical field. For computations it is better to call the parametric
!functions directly, they work in terms of hbar rather than x to
!prevent unnecessary numerical inversions
  function nmlfi3_x_star(xi,p,xend,w,lnRhoReh,Pstar,bfoldstar)
    implicit none
    real(kp) :: nmlfi3_x_star
    real(kp), intent(in) :: xi,p,xend,w,lnRhoReh,Pstar
    real(kp), intent(out), optional :: bfoldstar
    
    real(kp) :: hbarstar, hbarend

    hbarend = nmlfi_hbar(xend,xi)
    hbarstar = nmlfi3_hbar_star(xi,p,hbarend,w,lnRhoReh,Pstar,bfoldstar)

    nmlfi3_x_star = nmlfi_x(hbarstar,xi)
    
  end function nmlfi3_x_star


  function nmlfi3_x_rrad(xi,p,xend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: nmlfi3_x_rrad
    real(kp), intent(in) :: xi,p,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar
    
    real(kp) :: hbarstar,hbarend

    hbarend = nmlfi_hbar(xend,xi)
    hbarstar = nmlfi3_hbar_rrad(xi,p,hbarend,lnRrad,Pstar,bfoldstar)

    nmlfi3_x_rrad = nmlfi_x(hbarstar,xi)

  end function nmlfi3_x_rrad

    
  function nmlfi3_x_rreh(xi,p,xend,lnRreh,bfoldstar)
    implicit none
    real(kp) :: nmlfi3_x_rreh
    real(kp), intent(in) :: xi,p,xend,lnRreh
    real(kp), intent(out), optional :: bfoldstar
    
    real(kp) :: hbarstar,hbarend

    hbarend = nmlfi_hbar(xend,xi)
    hbarstar = nmlfi3_hbar_rreh(xi,p,hbarend,lnRreh,bfoldstar)

    nmlfi3_x_rreh = nmlfi_x(hbarstar,xi)
    
  end function nmlfi3_x_rreh
  
  
  
end module nmlficomreh

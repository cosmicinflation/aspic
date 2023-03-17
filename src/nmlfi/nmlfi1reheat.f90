!Reheating functions for Non-Minimal Large Field Inflation 1

module nmlfi1reheat
  use infprec, only :  kp
  use nmlficommon, only : nmlfi_x, nmlfi_hbar, hbarBig, hbarSmall
  use nmlficommon, only : nmlfi_hbar_potmax

  use nmlficomreh, only : nmlfi_hbar_rreh, nmlfi_hbar_star, nmlfi_hbar_rrad
  use nmlficomreh, only : nmlfi_lnrhoreh_max

  use nmlfi1sr, only : nmlfi1_check_params
  
  implicit none

  private

!take as input the parametric field hbar  
  public nmlfi1_parametric_lnrhoreh_max
  public nmlfi1_hbar_star, nmlfi1_hbar_rrad, nmlfi1_hbar_rreh


!for convenience only
  public nmlfi1_x_star, nmlfi1_x_rrad, nmlfi1_x_rreh, nmlfi1_lnrhoreh_max

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!the parametric reheating functions, the ones you should use !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function nmlfi1_hbar_star(xi,p,hbarend,w,lnRhoReh,Pstar,bfoldstar)
    implicit none
    real(kp) :: nmlfi1_hbar_star
    real(kp), intent(in) :: xi,p,hbarend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: hbarmin, hbarmax

    if (.not.nmlfi1_check_params(p,xi)) then
       write(*,*)'xi= p= ',xi,p
       stop 'nmlfi1_hbar_star: NMLFI1 does not exist for these parameters!'
    endif
    
    hbarmin = hbarend
    if (p.ge.4._kp) then
       hbarmax = hbarBig
    else
       hbarmax = nmlfi_hbar_potmax(xi,p)
    endif
    
    nmlfi1_hbar_star = nmlfi_hbar_star(xi,p,w,lnRhoReh,Pstar,hbarend,hbarmin,hbarmax,bfoldstar)
    
  end function nmlfi1_hbar_star



  function nmlfi1_hbar_rrad(xi,p,hbarend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: nmlfi1_hbar_rrad
    real(kp), intent(in) :: xi,p,hbarend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: hbarmin, hbarmax

    if (.not.nmlfi1_check_params(p,xi)) then
       write(*,*)'xi= p= ',xi,p
       stop 'nmlfi1_hbar_rrad: NMLFI1 does not exist for these parameters!'
    endif
    
    hbarmin = hbarend
    if (p.ge.4._kp) then
       hbarmax = hbarBig
    else
       hbarmax = nmlfi_hbar_potmax(xi,p)
    endif

    nmlfi1_hbar_rrad = nmlfi_hbar_rrad(xi,p,lnRrad,Pstar,hbarend,hbarmin,hbarmax,bfoldstar)

  end function nmlfi1_hbar_rrad
  


  function nmlfi1_hbar_rreh(xi,p,hbarend,lnRreh,bfoldstar)
    implicit none
    real(kp) :: nmlfi1_hbar_rreh
    real(kp), intent(in) :: xi,p,hbarend,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp) :: hbarmin, hbarmax

    if (.not.nmlfi1_check_params(p,xi)) then
       write(*,*)'xi= p= ',xi,p
       stop 'nmlfi1_hbar_rreh: NMLFI1 does not exist for these parameters!'
    endif
    
    hbarmin = hbarend
    if (p.ge.4._kp) then
       hbarmax = hbarBig
    else
       hbarmax = nmlfi_hbar_potmax(xi,p)
    endif

    nmlfi1_hbar_rreh = nmlfi_hbar_rreh(xi,p,lnRreh,hbarend,hbarmin,hbarmax,bfoldstar)
        
  end function nmlfi1_hbar_rreh

  
    
!jordan frame
  function nmlfi1_parametric_lnrhoreh_max(xi,p,hbarend,Pstar)
    implicit none
    real(kp) :: nmlfi1_parametric_lnrhoreh_max
    real(kp), intent(in) :: xi,p,hbarend,Pstar
    real(kp) :: hbarmin, hbarmax
        
    hbarmin = hbarend
    if (p.ge.4._kp) then
       hbarmax = hbarBig
    else
       hbarmax = nmlfi_hbar_potmax(xi,p)
    endif
    
    nmlfi1_parametric_lnrhoreh_max = nmlfi_lnrhoreh_max(xi,p,Pstar,hbarend,hbarmin,hbarmax)

  end function nmlfi1_parametric_lnrhoreh_max




!Convenience functions: inefficient wrapper taking as input the
!canonical field. For computations it is better to call the parametric
!functions directly, they work in terms of hbar rather than x to
!prevent unnecessary numerical inversions
  function nmlfi1_x_star(xi,p,xend,w,lnRhoReh,Pstar,bfoldstar)
    implicit none
    real(kp) :: nmlfi1_x_star
    real(kp), intent(in) :: xi,p,xend,w,lnRhoReh,Pstar
    real(kp), intent(out), optional :: bfoldstar
    
    real(kp) :: hbarstar, hbarend

    hbarend = nmlfi_hbar(xend,xi)
    hbarstar = nmlfi1_hbar_star(xi,p,hbarend,w,lnRhoReh,Pstar,bfoldstar)

    nmlfi1_x_star = nmlfi_x(hbarstar,xi)
    
  end function nmlfi1_x_star


  function nmlfi1_x_rrad(xi,p,xend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: nmlfi1_x_rrad
    real(kp), intent(in) :: xi,p,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar
    
    real(kp) :: hbarstar,hbarend

    hbarend = nmlfi_hbar(xend,xi)
    hbarstar = nmlfi1_hbar_rrad(xi,p,hbarend,lnRrad,Pstar,bfoldstar)

    nmlfi1_x_rrad = nmlfi_x(hbarstar,xi)

  end function nmlfi1_x_rrad

    
  function nmlfi1_x_rreh(xi,p,xend,lnRreh,bfoldstar)
    implicit none
    real(kp) :: nmlfi1_x_rreh
    real(kp), intent(in) :: xi,p,xend,lnRreh
    real(kp), intent(out), optional :: bfoldstar
    
    real(kp) :: hbarstar,hbarend

    hbarend = nmlfi_hbar(xend,xi)
    hbarstar = nmlfi1_hbar_rreh(xi,p,hbarend,lnRreh,bfoldstar)

    nmlfi1_x_rreh = nmlfi_x(hbarstar,xi)
    
  end function nmlfi1_x_rreh


  function nmlfi1_lnrhoreh_max(xi,p,xend,Pstar)
    implicit none
    real(kp) :: nmlfi1_lnrhoreh_max
    real(kp), intent(in) :: xi,p,xend,Pstar
    real(kp) :: hbarmin, hbarmax, hbarend
    
    hbarend = nmlfi_hbar(xend,xi)
  
    nmlfi1_lnrhoreh_max = nmlfi1_parametric_lnrhoreh_max(xi,p,hbarend,Pstar)

  end function nmlfi1_lnrhoreh_max

  
  
  
end module nmlfi1reheat

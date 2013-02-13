!radiatively corrected massive inflation reheating functions in the slow-roll

module rcmireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use rcmisr, only : rcmi_epsilon_one, rcmi_epsilon_two
  use rcmisr, only : rcmi_norm_potential, rcmi_x_potmax
  use rcmisr, only : rcmi_x_endinf, rcmi_efold_primitive
  use specialinf, only : lambert
  implicit none

  private

  public rcmi_x_star, rcmi_lnrhoend 

contains

!returns x given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the correspoding bfoldstar
  function rcmi_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: rcmi_x_star
    real(kp), intent(in) :: alpha,lnRhoReh,w,Pstar
    real(kp), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x,xPotMax
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: rcmiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = rcmi_x_endinf(alpha)
    xPotMax = rcmi_x_potmax(alpha)
    epsOneEnd = rcmi_epsilon_one(xEnd,alpha)
    potEnd = rcmi_norm_potential(xEnd,alpha)
    primEnd = rcmi_efold_primitive(xEnd,alpha)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    rcmiData%real1 = alpha
    rcmiData%real2 = w
    rcmiData%real3 = calF + primEnd

    mini = xEnd
    maxi = min(xPotMax,1._kp/epsilon(1._kp))

    x = zbrent(find_rcmi_x_star,mini,maxi,tolzbrent,rcmiData)
    rcmi_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (rcmi_efold_primitive(x,alpha) - primEnd)
    endif
    
  end function rcmi_x_star

  function find_rcmi_x_star(x,rcmiData)   
    implicit none
    real(kp) :: find_rcmi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rcmiData

    real(kp) :: primStar,alpha,w,CalFplusprimEnd,potStar,epsOneStar

    alpha = rcmiData%real1
    w = rcmiData%real2
    CalFplusprimEnd = rcmiData%real3

    primStar = rcmi_efold_primitive(x,alpha)
    epsOneStar = rcmi_epsilon_one(x,alpha)
    potStar = rcmi_norm_potential(x,alpha)

    find_rcmi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)

  end function find_rcmi_x_star



  function rcmi_lnrhoend(alpha,Pstar) 
    implicit none
    real(kp) :: rcmi_lnrhoend
    real(kp), intent(in) :: alpha,Pstar

    real(kp) :: xEnd, potEnd, epsEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp), parameter :: wrad = 1._kp/3._kp
    real(kp), parameter :: junk = 0._kp

    real(kp) :: lnRhoEnd
    
    xEnd = rcmi_x_endinf(alpha)      
    potEnd  = rcmi_norm_potential(xEnd,alpha)
    epsEnd = rcmi_epsilon_one(xEnd,alpha)

!   Trick to return x such that rho_reh=rho_end
       
    x = rcmi_x_star(alpha,wrad,junk,Pstar)    
    potStar = rcmi_norm_potential(x,alpha)
    epsOneStar = rcmi_epsilon_one(x,alpha)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'rcmi_lnrhoend: slow-roll violated!'

    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsEnd,potEnd/potStar)

    rcmi_lnrhoend = lnRhoEnd

  end function rcmi_lnrhoend

  

    
end module rcmireheat

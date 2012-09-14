!Valley Hybrid inflation reheating functions in the slow-roll approximations

module vhireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use vhisr, only : vhi_epsilon_one, vhi_epsilon_two, vhi_epsilon_three
  use vhisr, only : vhi_norm_potential,vhi_efold_primitive,vhi_xend_max,vhi_xend_min
  implicit none

  private

  public vhi_x_star, vhi_lnrhoend 

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function vhi_x_star(p,mu,xEnd,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: vhi_x_star
    real(kp), intent(in) :: p,mu,xEnd,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: vhiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    

    epsOneEnd = vhi_epsilon_one(xEnd,p,mu)
    potEnd = vhi_norm_potential(xEnd,p)

    primEnd = vhi_efold_primitive(xEnd,p,mu)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    vhiData%real1 = p
    vhiData%real2 = mu
    vhiData%real3 = w
    vhiData%real4 = calF + primEnd

    mini = xEnd*(1._kp+epsilon(1._kp))
    maxi= vhi_xend_max(p,mu)*(1._kp-epsilon(1._kp))

    x = zbrent(find_vhi_x_star,mini,maxi,tolzbrent,vhiData)
    vhi_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (vhi_efold_primitive(x,p,mu) - primEnd)
    endif

  end function vhi_x_star

  function find_vhi_x_star(x,vhiData)   
    implicit none
    real(kp) :: find_vhi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: vhiData

    real(kp) :: primStar,p,mu,w,CalFplusprimEnd,potStar,epsOneStar

    p=vhiData%real1
    mu=vhiData%real2
    w = vhiData%real3
    CalFplusprimEnd = vhiData%real4

    primStar = vhi_efold_primitive(x,p,mu)
    epsOneStar = vhi_epsilon_one(x,p,mu)
    potStar = vhi_norm_potential(x,p)

    find_vhi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_vhi_x_star



  function vhi_lnrhoend(p,mu,xEnd,Pstar) 
    implicit none
    real(kp) :: vhi_lnrhoend
    real(kp), intent(in) :: p,mu,xEnd,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd


    potEnd  = vhi_norm_potential(xEnd,p)

    epsOneEnd = vhi_epsilon_one(xEnd,p,mu)

!   Trick to return x such that rho_reh=rho_end

    x = vhi_x_star(p,mu,xEnd,wrad,junk,Pstar)  

 
    potStar = vhi_norm_potential(x,p)
    epsOneStar = vhi_epsilon_one(x,p,mu)

    
    if (.not.slowroll_validity(epsOneStar)) stop 'vhi_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    vhi_lnrhoend = lnRhoEnd

  end function vhi_lnrhoend

  
end module vhireheat

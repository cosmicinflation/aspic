!sneutrino supersymmetric 6 reheating functions in the slow-roll approximations

module ssi6reheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use ssi6sr, only : ssi6_epsilon_one, ssi6_epsilon_two, ssi6_epsilon_three
  use ssi6sr, only : ssi6_norm_potential
  use ssi6sr, only : ssi6_x_endinf, ssi6_efold_primitive
  use cosmopar, only : QrmsOverT

  implicit none

  private

  public ssi6_x_star, ssi6_lnrhoend

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function ssi6_x_star(alpha,beta,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: ssi6_x_star
    real(kp), intent(in) :: alpha,beta,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd
    type(transfert) :: ssi6Data

  
    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif

    xEnd=ssi6_x_endinf(alpha,beta)
    epsOneEnd = ssi6_epsilon_one(xEnd,alpha,beta)
    potEnd = ssi6_norm_potential(xEnd,alpha,beta)

    primEnd = ssi6_efold_primitive(xEnd,alpha,beta)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    ssi6Data%real1 = alpha
    ssi6Data%real2 = beta
    ssi6Data%real3 = w
    ssi6Data%real4 = calF + primEnd

    mini = ssi6_x_endinf(alpha,beta)*(1._kp+epsilon(1._kp))
    maxi = 10._kp**(6._kp)*mini

    x = zbrent(find_ssi6_x_star,mini,maxi,tolzbrent,ssi6Data)
    ssi6_x_star = x

!   print*,'ssi6_x_star:  xEnd=',xEnd,'  potEnd=',potEnd,'   epsOneEnd=',epsOneEnd, &
!       '   primEnd=',primEnd,'   xstar=',x
!    pause

    if (present(bfoldstar)) then
       bfoldstar = - (ssi6_efold_primitive(x,alpha,beta) - primEnd)
    endif

  end function ssi6_x_star

  function find_ssi6_x_star(x,ssi6Data)   
    implicit none
    real(kp) :: find_ssi6_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: ssi6Data

    real(kp) :: primStar,alpha,beta,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=ssi6Data%real1
    beta=ssi6Data%real2
    w = ssi6Data%real3
    CalFplusprimEnd = ssi6Data%real4

    primStar = ssi6_efold_primitive(x,alpha,beta)
    epsOneStar = ssi6_epsilon_one(x,alpha,beta)
    potStar = ssi6_norm_potential(x,alpha,beta)

    find_ssi6_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_ssi6_x_star



  function ssi6_lnrhoend(alpha,beta,Pstar) 
    implicit none
    real(kp) :: ssi6_lnrhoend
    real(kp), intent(in) :: alpha,beta,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = ssi6_x_endinf(alpha,beta)

    potEnd  = ssi6_norm_potential(xEnd,alpha,beta)

    epsOneEnd = ssi6_epsilon_one(xEnd,alpha,beta)

!   print*,'ssi6_lnrhoend:  xEnd=',xEnd,'  potEnd=',potEnd,'   epsOneEnd=',epsOneEnd
!    pause


!   Trick to return x such that rho_reh=rho_end

    x = ssi6_x_star(alpha,beta,wrad,junk,Pstar)  


    potStar = ssi6_norm_potential(x,alpha,beta)
    epsOneStar = ssi6_epsilon_one(x,alpha,beta)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'ssi6_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    ssi6_lnrhoend = lnRhoEnd

  end function ssi6_lnrhoend

  
end module ssi6reheat

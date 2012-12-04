!sneutrino supersymmetric 3 reheating functions in the slow-roll approximations

module ssi3reheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use ssi3sr, only : ssi3_epsilon_one, ssi3_epsilon_two, ssi3_epsilon_three
  use ssi3sr, only : ssi3_norm_potential
  use ssi3sr, only : ssi3_x_endinf, ssi3_efold_primitive
  use ssicommon, only : ssi3456_x_Vprime_Equals_0
  use cosmopar, only : QrmsOverT

  implicit none

  private

  public ssi3_x_star, ssi3_lnrhoend

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function ssi3_x_star(alpha,beta,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: ssi3_x_star
    real(kp), intent(in) :: alpha,beta,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd
    type(transfert) :: ssi3Data

  
    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif

    xEnd=ssi3_x_endinf(alpha,beta)
    epsOneEnd = ssi3_epsilon_one(xEnd,alpha,beta)
    potEnd = ssi3_norm_potential(xEnd,alpha,beta)

    primEnd = ssi3_efold_primitive(xEnd,alpha,beta)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    ssi3Data%real1 = alpha
    ssi3Data%real2 = beta
    ssi3Data%real3 = w
    ssi3Data%real4 = calF + primEnd

    mini = ssi3_x_endinf(alpha,beta)*(1._kp+epsilon(1._kp))
    maxi = ssi3456_x_Vprime_Equals_0(alpha,beta)*(1._kp-epsilon(1._kp))

    x = zbrent(find_ssi3_x_star,mini,maxi,tolzbrent,ssi3Data)
    ssi3_x_star = x

!   print*,'ssi3_x_star:  xEnd=',xEnd,'  potEnd=',potEnd,'   epsOneEnd=',epsOneEnd, &
!       '   primEnd=',primEnd,'   xstar=',x
!    pause

    if (present(bfoldstar)) then
       bfoldstar = - (ssi3_efold_primitive(x,alpha,beta) - primEnd)
    endif

  end function ssi3_x_star

  function find_ssi3_x_star(x,ssi3Data)   
    implicit none
    real(kp) :: find_ssi3_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: ssi3Data

    real(kp) :: primStar,alpha,beta,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=ssi3Data%real1
    beta=ssi3Data%real2
    w = ssi3Data%real3
    CalFplusprimEnd = ssi3Data%real4

    primStar = ssi3_efold_primitive(x,alpha,beta)
    epsOneStar = ssi3_epsilon_one(x,alpha,beta)
    potStar = ssi3_norm_potential(x,alpha,beta)

    find_ssi3_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_ssi3_x_star



  function ssi3_lnrhoend(alpha,beta,Pstar) 
    implicit none
    real(kp) :: ssi3_lnrhoend
    real(kp), intent(in) :: alpha,beta,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = ssi3_x_endinf(alpha,beta)

    potEnd  = ssi3_norm_potential(xEnd,alpha,beta)

    epsOneEnd = ssi3_epsilon_one(xEnd,alpha,beta)

!   print*,'ssi3_lnrhoend:  xEnd=',xEnd,'  potEnd=',potEnd,'   epsOneEnd=',epsOneEnd
!    pause


!   Trick to return x such that rho_reh=rho_end

    x = ssi3_x_star(alpha,beta,wrad,junk,Pstar)  


    potStar = ssi3_norm_potential(x,alpha,beta)
    epsOneStar = ssi3_epsilon_one(x,alpha,beta)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'ssi3_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    ssi3_lnrhoend = lnRhoEnd

  end function ssi3_lnrhoend

  
end module ssi3reheat

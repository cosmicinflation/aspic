!sneutrino supersymmetric 2 reheating functions in the slow-roll approximations

module ssi2reheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use ssi2sr, only : ssi2_epsilon_one, ssi2_epsilon_two, ssi2_epsilon_three
  use ssi2sr, only : ssi2_norm_potential
  use ssi2sr, only : ssi2_x_endinf, ssi2_efold_primitive
  use cosmopar, only : QrmsOverT
  implicit none

  private

  public ssi2_x_star, ssi2_lnrhoend

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function ssi2_x_star(alpha,beta,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: ssi2_x_star
    real(kp), intent(in) :: alpha,beta,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd
    type(transfert) :: ssi2Data

  
    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif

    xEnd=ssi2_x_endinf(alpha,beta)
    epsOneEnd = ssi2_epsilon_one(xEnd,alpha,beta)
    potEnd = ssi2_norm_potential(xEnd,alpha,beta)

    primEnd = ssi2_efold_primitive(xEnd,alpha,beta)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    ssi2Data%real1 = alpha
    ssi2Data%real2 = beta
    ssi2Data%real3 = w
    ssi2Data%real4 = calF + primEnd

    mini = epsilon(1._kp)
    maxi = ssi2_x_endinf(alpha,beta)*(1._kp-epsilon(1._kp))

    x = zbrent(find_ssi2_x_star,mini,maxi,tolzbrent,ssi2Data)
    ssi2_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (ssi2_efold_primitive(x,alpha,beta) - primEnd)
    endif

  end function ssi2_x_star

  function find_ssi2_x_star(x,ssi2Data)   
    implicit none
    real(kp) :: find_ssi2_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: ssi2Data

    real(kp) :: primStar,alpha,beta,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=ssi2Data%real1
    beta=ssi2Data%real2
    w = ssi2Data%real3
    CalFplusprimEnd = ssi2Data%real4

    primStar = ssi2_efold_primitive(x,alpha,beta)
    epsOneStar = ssi2_epsilon_one(x,alpha,beta)
    potStar = ssi2_norm_potential(x,alpha,beta)

    find_ssi2_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_ssi2_x_star



  function ssi2_lnrhoend(alpha,beta,Pstar) 
    implicit none
    real(kp) :: ssi2_lnrhoend
    real(kp), intent(in) :: alpha,beta,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = ssi2_x_endinf(alpha,beta)

    potEnd  = ssi2_norm_potential(xEnd,alpha,beta)

    epsOneEnd = ssi2_epsilon_one(xEnd,alpha,beta)

!    print*,'ssi2_lnrhoend:  xEnd=',xEnd,'  potEnd=',potEnd,'   epsOneEnd=',epsOneEnd
!    pause


!   Trick to return x such that rho_reh=rho_end

    x = ssi2_x_star(alpha,beta,wrad,junk,Pstar)  


    potStar = ssi2_norm_potential(x,alpha,beta)
    epsOneStar = ssi2_epsilon_one(x,alpha,beta)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'ssi2_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    ssi2_lnrhoend = lnRhoEnd

  end function ssi2_lnrhoend

  
end module ssi2reheat

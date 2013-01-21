!sneutrino supersymmetric 4 reheating functions in the slow-roll approximations

module ssi4reheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use ssi4sr, only : ssi4_epsilon_one, ssi4_epsilon_two, ssi4_epsilon_three
  use ssi4sr, only : ssi4_norm_potential, ssi4_x_potmax
  use ssi4sr, only : ssi4_x_endinf, ssi4_efold_primitive
  use cosmopar, only : QrmsOverT
  implicit none

  private

  public ssi4_x_star, ssi4_lnrhoend

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function ssi4_x_star(alpha,beta,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: ssi4_x_star
    real(kp), intent(in) :: alpha,beta,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd
    type(transfert) :: ssi4Data

  
    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif

    xEnd=ssi4_x_endinf(alpha,beta)
    epsOneEnd = ssi4_epsilon_one(xEnd,alpha,beta)
    potEnd = ssi4_norm_potential(xEnd,alpha,beta)

    primEnd = ssi4_efold_primitive(xEnd,alpha,beta)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    ssi4Data%real1 = alpha
    ssi4Data%real2 = beta
    ssi4Data%real3 = w
    ssi4Data%real4 = calF + primEnd


    mini = ssi4_x_potmax(alpha,beta)*(1._kp+epsilon(1._kp))
    maxi = ssi4_x_endinf(alpha,beta)*(1._kp-epsilon(1._kp))

    x = zbrent(find_ssi4_x_star,mini,maxi,tolzbrent,ssi4Data)
    ssi4_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (ssi4_efold_primitive(x,alpha,beta) - primEnd)
    endif

  end function ssi4_x_star

  function find_ssi4_x_star(x,ssi4Data)   
    implicit none
    real(kp) :: find_ssi4_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: ssi4Data

    real(kp) :: primStar,alpha,beta,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=ssi4Data%real1
    beta=ssi4Data%real2
    w = ssi4Data%real3
    CalFplusprimEnd = ssi4Data%real4

    primStar = ssi4_efold_primitive(x,alpha,beta)
    epsOneStar = ssi4_epsilon_one(x,alpha,beta)
    potStar = ssi4_norm_potential(x,alpha,beta)

    find_ssi4_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_ssi4_x_star



  function ssi4_lnrhoend(alpha,beta,Pstar) 
    implicit none
    real(kp) :: ssi4_lnrhoend
    real(kp), intent(in) :: alpha,beta,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = ssi4_x_endinf(alpha,beta)

    potEnd  = ssi4_norm_potential(xEnd,alpha,beta)

    epsOneEnd = ssi4_epsilon_one(xEnd,alpha,beta)


!   Trick to return x such that rho_reh=rho_end

    x = ssi4_x_star(alpha,beta,wrad,junk,Pstar)  


    potStar = ssi4_norm_potential(x,alpha,beta)
    epsOneStar = ssi4_epsilon_one(x,alpha,beta)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'ssi4_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    ssi4_lnrhoend = lnRhoEnd

  end function ssi4_lnrhoend

  
end module ssi4reheat

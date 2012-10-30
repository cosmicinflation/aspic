!sneutrino supersymmetric 1 reheating functions in the slow-roll approximations

module ssi1reheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use ssi1sr, only : ssi1_epsilon_one, ssi1_epsilon_two, ssi1_epsilon_three
  use ssi1sr, only : ssi1_norm_potential
  use ssi1sr, only : ssi1_x_endinf, ssi1_efold_primitive
  use cosmopar, only : QrmsOverT
  implicit none

  private

  public ssi1_x_star, ssi1_lnrhoend

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function ssi1_x_star(alpha,beta,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: ssi1_x_star
    real(kp), intent(in) :: alpha,beta,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd
    type(transfert) :: ssi1Data

  
    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif

    xEnd=ssi1_x_endinf(alpha,beta)
    epsOneEnd = ssi1_epsilon_one(xEnd,alpha,beta)
    potEnd = ssi1_norm_potential(xEnd,alpha,beta)

    primEnd = ssi1_efold_primitive(xEnd,alpha,beta)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    ssi1Data%real1 = alpha
    ssi1Data%real2 = beta
    ssi1Data%real3 = w
    ssi1Data%real4 = calF + primEnd

    mini = xEnd*(1._kp+epsilon(1._kp))
    maxi = mini/epsilon(1._kp)

    x = zbrent(find_ssi1_x_star,mini,maxi,tolzbrent,ssi1Data)
    ssi1_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (ssi1_efold_primitive(x,alpha,beta) - primEnd)
    endif

  end function ssi1_x_star

  function find_ssi1_x_star(x,ssi1Data)   
    implicit none
    real(kp) :: find_ssi1_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: ssi1Data

    real(kp) :: primStar,alpha,beta,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=ssi1Data%real1
    beta=ssi1Data%real2
    w = ssi1Data%real3
    CalFplusprimEnd = ssi1Data%real4

    primStar = ssi1_efold_primitive(x,alpha,beta)
    epsOneStar = ssi1_epsilon_one(x,alpha,beta)
    potStar = ssi1_norm_potential(x,alpha,beta)

    find_ssi1_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_ssi1_x_star



  function ssi1_lnrhoend(alpha,beta,Pstar) 
    implicit none
    real(kp) :: ssi1_lnrhoend
    real(kp), intent(in) :: alpha,beta,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = ssi1_x_endinf(alpha,beta)

    potEnd  = ssi1_norm_potential(xEnd,alpha,beta)

    epsOneEnd = ssi1_epsilon_one(xEnd,alpha,beta)


!   Trick to return x such that rho_reh=rho_end

    x = ssi1_x_star(alpha,beta,wrad,junk,Pstar)  


    potStar = ssi1_norm_potential(x,alpha,beta)
    epsOneStar = ssi1_epsilon_one(x,alpha,beta)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'ssi1_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    ssi1_lnrhoend = lnRhoEnd

  end function ssi1_lnrhoend

  
end module ssi1reheat

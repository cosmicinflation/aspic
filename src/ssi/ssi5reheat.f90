!sneutrino supersymmetric 5 reheating functions in the slow-roll approximations

module ssi5reheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use ssi5sr, only : ssi5_epsilon_one, ssi5_epsilon_two, ssi5_epsilon_three
  use ssi5sr, only : ssi5_norm_potential
  use ssi5sr, only : ssi5_x_endinf, ssi5_efold_primitive
  use cosmopar, only : QrmsOverT

  implicit none

  private

  public ssi5_x_star, ssi5_lnrhoend

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function ssi5_x_star(alpha,beta,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: ssi5_x_star
    real(kp), intent(in) :: alpha,beta,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd
    type(transfert) :: ssi5Data

  
    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif

    xEnd=ssi5_x_endinf(alpha,beta)
    epsOneEnd = ssi5_epsilon_one(xEnd,alpha,beta)
    potEnd = ssi5_norm_potential(xEnd,alpha,beta)

    primEnd = ssi5_efold_primitive(xEnd,alpha,beta)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    ssi5Data%real1 = alpha
    ssi5Data%real2 = beta
    ssi5Data%real3 = w
    ssi5Data%real4 = calF + primEnd

    mini = epsilon(1._kp)
    maxi = ssi5_x_endinf(alpha,beta)*(1._kp-epsilon(1._kp))

    x = zbrent(find_ssi5_x_star,mini,maxi,tolzbrent,ssi5Data)
    ssi5_x_star = x

!   print*,'ssi5_x_star:  xEnd=',xEnd,'  potEnd=',potEnd,'   epsOneEnd=',epsOneEnd, &
!       '   primEnd=',primEnd,'   xstar=',x
!    pause

    if (present(bfoldstar)) then
       bfoldstar = - (ssi5_efold_primitive(x,alpha,beta) - primEnd)
    endif

  end function ssi5_x_star

  function find_ssi5_x_star(x,ssi5Data)   
    implicit none
    real(kp) :: find_ssi5_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: ssi5Data

    real(kp) :: primStar,alpha,beta,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=ssi5Data%real1
    beta=ssi5Data%real2
    w = ssi5Data%real3
    CalFplusprimEnd = ssi5Data%real4

    primStar = ssi5_efold_primitive(x,alpha,beta)
    epsOneStar = ssi5_epsilon_one(x,alpha,beta)
    potStar = ssi5_norm_potential(x,alpha,beta)

    find_ssi5_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_ssi5_x_star



  function ssi5_lnrhoend(alpha,beta,Pstar) 
    implicit none
    real(kp) :: ssi5_lnrhoend
    real(kp), intent(in) :: alpha,beta,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = ssi5_x_endinf(alpha,beta)

    potEnd  = ssi5_norm_potential(xEnd,alpha,beta)

    epsOneEnd = ssi5_epsilon_one(xEnd,alpha,beta)

!   print*,'ssi5_lnrhoend:  xEnd=',xEnd,'  potEnd=',potEnd,'   epsOneEnd=',epsOneEnd
!    pause


!   Trick to return x such that rho_reh=rho_end

    x = ssi5_x_star(alpha,beta,wrad,junk,Pstar)  


    potStar = ssi5_norm_potential(x,alpha,beta)
    epsOneStar = ssi5_epsilon_one(x,alpha,beta)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'ssi5_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    ssi5_lnrhoend = lnRhoEnd

  end function ssi5_lnrhoend

  
end module ssi5reheat

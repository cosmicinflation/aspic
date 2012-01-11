!intermediate inflation reheating functions in the slow-roll approximations

module iireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use iisr, only : ii_epsilon_one, ii_epsilon_two, ii_epsilon_three
  use iisr, only : ii_norm_potential
  use iisr, only : ii_x_endinf, ii_efold_primitive
  implicit none

  private

  public ii_x_star, ii_lnrhoend 

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function ii_x_star(beta,xendinf,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: ii_x_star
    real(kp), intent(in) :: beta,xendinf,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: iiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = ii_x_endinf(beta,xendinf)
    epsOneEnd = ii_epsilon_one(xEnd,beta)
    potEnd = ii_norm_potential(xEnd,beta)
    primEnd = ii_efold_primitive(xEnd,beta)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    iiData%real1 = beta 
    iiData%real2 = xEnd
    iiData%real3 = w
    iiData%real4 = calF + primEnd

    mini = sqrt(beta/2._kp*(1._kp+beta/3._kp+sqrt(1._kp+4._kp*beta/9._kp)))
    maxi = xEnd


    x = zbrent(find_ii_x_star,mini,maxi,tolzbrent,iiData)
    ii_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (ii_efold_primitive(x,beta) - primEnd)
    endif

  end function ii_x_star

  function find_ii_x_star(x,iiData)   
    implicit none
    real(kp) :: find_ii_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: iiData

    real(kp) :: primStar,beta,xEnd,w,CalFplusprimEnd,potStar,epsOneStar

    beta=iiData%real1
    xEnd=iiData%real2
    w = iiData%real3
    CalFplusprimEnd = iiData%real4

    primStar = ii_efold_primitive(x,beta)
    epsOneStar = ii_epsilon_one(x,beta)
    potStar = ii_norm_potential(x,beta)

    find_ii_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_ii_x_star



  function ii_lnrhoend(beta,xendinf,Pstar) 
    implicit none
    real(kp) :: ii_lnrhoend
    real(kp), intent(in) :: beta,xendinf,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = ii_x_endinf(beta,xendinf)
    potEnd  = ii_norm_potential(xEnd,beta)
    epsOneEnd = ii_epsilon_one(xEnd,beta)


!   Trick to return x such that rho_reh=rho_end

    x = ii_x_star(beta,xEnd,wrad,junk,Pstar)    
    potStar = ii_norm_potential(x,beta)
    epsOneStar = ii_epsilon_one(x,beta)


    
    if (.not.slowroll_validity(epsOneStar)) stop 'ii_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    ii_lnrhoend = lnRhoEnd

  end function ii_lnrhoend

  
end module iireheat

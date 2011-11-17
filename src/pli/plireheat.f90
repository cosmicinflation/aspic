!power law inflation reheating functions in the slow-roll approximations

module plireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use plisr, only : pli_epsilon_one, pli_epsilon_two, pli_epsilon_three
  use plisr, only : pli_norm_potential
  use plisr, only : pli_efold_primitive
  implicit none

  private

  public pli_x_star, pli_lnrhoend 

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function pli_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: pli_x_star
    real(kp), intent(in) :: alpha,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
 !   real(kp), parameter :: junk_xend=1000._kp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: pliData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
!    xEnd = junk_xend
     xEnd = 100.
    epsOneEnd = pli_epsilon_one(xEnd,alpha)
    potEnd = pli_norm_potential(xEnd,alpha)
    primEnd = pli_efold_primitive(xEnd,alpha)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    pliData%real1 = alpha    
    pliData%real2 = w
    pliData%real3 = calF + primEnd

!    mini = epsilon(1._kp)
    mini = 1.
    maxi = xend


    x = zbrent(find_pli_x_star,mini,maxi,tolzbrent,pliData)
    pli_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (pli_efold_primitive(x,alpha) - primEnd)
    endif

  end function pli_x_star

  function find_pli_x_star(x,pliData)   
    implicit none
    real(kp) :: find_pli_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: pliData

    real(kp) :: primStar,alpha,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=pliData%real1
    w = pliData%real2
    CalFplusprimEnd = pliData%real3

    primStar = pli_efold_primitive(x,alpha)
    epsOneStar = pli_epsilon_one(x,alpha)
    potStar = pli_norm_potential(x,alpha)

    find_pli_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_pli_x_star



  function pli_lnrhoend(alpha,Pstar) 
    implicit none
    real(kp) :: pli_lnrhoend
    real(kp), intent(in) :: alpha,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp
!    real(kp),parameter :: junk_xend=100_kp

    real(kp) :: lnRhoEnd
    
    !xEnd = junk_xend
    xEnd=100._kp
    potEnd  = pli_norm_potential(xEnd,alpha)
    epsOneEnd = pli_epsilon_one(xEnd,alpha)

!   Trick to return x such that rho_reh=rho_end

    x = pli_x_star(alpha,wrad,junk,Pstar)    
    potStar = pli_norm_potential(x,alpha)
    epsOneStar = pli_epsilon_one(x,alpha)

    
    if (.not.slowroll_validity(epsOneStar)) stop 'pli_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    pli_lnrhoend = lnRhoEnd

  end function pli_lnrhoend

  
end module plireheat

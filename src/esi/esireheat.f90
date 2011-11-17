!exponential SUSY inflation reheating functions in the slow-roll approximations

module esireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use esisr, only : esi_epsilon_one, esi_epsilon_two, esi_epsilon_three
  use esisr, only : esi_norm_potential
  use esisr, only : esi_x_endinf, esi_efold_primitive
  implicit none

  private

  public esi_x_star, esi_lnrhoend 

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function esi_x_star(q,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: esi_x_star
    real(kp), intent(in) :: q,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: esiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = esi_x_endinf(q)
    epsOneEnd = esi_epsilon_one(xEnd,q)
    potEnd = esi_norm_potential(xEnd,q)
    primEnd = esi_efold_primitive(xEnd,q)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    esiData%real1 = q    
    esiData%real2 = w
    esiData%real3 = calF + primEnd

    mini = xEnd
    maxi = 1._kp/epsilon(1._kp)

    x = zbrent(find_esi_x_star,mini,maxi,tolzbrent,esiData)
    esi_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (esi_efold_primitive(x,q) - primEnd)
    endif

  end function esi_x_star

  function find_esi_x_star(x,esiData)   
    implicit none
    real(kp) :: find_esi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: esiData

    real(kp) :: primStar,q,w,CalFplusprimEnd,potStar,epsOneStar

    q=esiData%real1
    w = esiData%real2
    CalFplusprimEnd = esiData%real3

    primStar = esi_efold_primitive(x,q)
    epsOneStar = esi_epsilon_one(x,q)
    potStar = esi_norm_potential(x,q)

    find_esi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_esi_x_star



  function esi_lnrhoend(q,Pstar) 
    implicit none
    real(kp) :: esi_lnrhoend
    real(kp), intent(in) :: q,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = esi_x_endinf(q)
    potEnd  = esi_norm_potential(xEnd,q)
    epsOneEnd = esi_epsilon_one(xEnd,q)

!   Trick to return x such that rho_reh=rho_end

    x = esi_x_star(q,wrad,junk,Pstar)    
    potStar = esi_norm_potential(x,q)
    epsOneStar = esi_epsilon_one(x,q)

    
    if (.not.slowroll_validity(epsOneStar)) stop 'esi_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    esi_lnrhoend = lnRhoEnd

  end function esi_lnrhoend

  
end module esireheat

!inverse monomial model reheating functions in the slow-roll approximations

module imireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use imisr, only : imi_epsilon_one, imi_epsilon_two, imi_epsilon_three
  use imisr, only : imi_norm_potential,imi_x_epsOne_Equal_One
  use imisr, only : imi_efold_primitive
  implicit none

  private

  public imi_x_star, imi_lnrhoend 

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function imi_x_star(p,xend,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: imi_x_star
    real(kp), intent(in) :: p,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: imiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = imi_epsilon_one(xEnd,p)
    potEnd = imi_norm_potential(xEnd,p)
    primEnd = imi_efold_primitive(xEnd,p)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    imiData%real1 = p    
    imiData%real2 = w
    imiData%real3 = calF + primEnd

    mini = imi_x_epsOne_Equal_One(p)
    maxi = xEnd

    x = zbrent(find_imi_x_star,mini,maxi,tolzbrent,imiData)
    imi_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (imi_efold_primitive(x,p) - primEnd)
    endif

  end function imi_x_star

  function find_imi_x_star(x,imiData)   
    implicit none
    real(kp) :: find_imi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: imiData

    real(kp) :: primStar,p,w,CalFplusprimEnd,potStar,epsOneStar

    p=imiData%real1
    w = imiData%real2
    CalFplusprimEnd = imiData%real3

    primStar = imi_efold_primitive(x,p)
    epsOneStar = imi_epsilon_one(x,p)
    potStar = imi_norm_potential(x,p)

    find_imi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_imi_x_star



  function imi_lnrhoend(p,xend,Pstar) 
    implicit none
    real(kp) :: imi_lnrhoend
    real(kp), intent(in) :: p,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    potEnd  = imi_norm_potential(xEnd,p)
    epsOneEnd = imi_epsilon_one(xEnd,p)

!   Trick to return x such that rho_reh=rho_end

    x = imi_x_star(p,xend,wrad,junk,Pstar)    
    potStar = imi_norm_potential(x,p)
    epsOneStar = imi_epsilon_one(x,p)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'imi_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    imi_lnrhoend = lnRhoEnd

  end function imi_lnrhoend

  
end module imireheat

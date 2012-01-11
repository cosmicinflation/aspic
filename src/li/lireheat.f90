!loop inflation reheating functions in the slow-roll approximations

module lireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use lisr, only : li_epsilon_one, li_epsilon_two, li_epsilon_three
  use lisr, only : li_norm_potential
  use lisr, only : li_x_endinf, li_efold_primitive
  implicit none

  private

  public li_x_star, li_lnrhoend 

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function li_x_star(alpha,mu,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: li_x_star
    real(kp), intent(in) :: alpha,mu,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: liData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = li_x_endinf(alpha,mu)
    epsOneEnd = li_epsilon_one(xEnd,alpha,mu)
    potEnd = li_norm_potential(xEnd,alpha,mu)
    primEnd = li_efold_primitive(xEnd,alpha,mu)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    liData%real1 = alpha
    liData%real2=mu   
    liData%real3 = w
    liData%real4 = calF + primEnd

    mini = xEnd
    maxi = xEnd*1000._kp

    x = zbrent(find_li_x_star,mini,maxi,tolzbrent,liData)
    li_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (li_efold_primitive(x,alpha,mu) - primEnd)
    endif

  end function li_x_star

  function find_li_x_star(x,liData)   
    implicit none
    real(kp) :: find_li_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: liData

    real(kp) :: primStar,alpha,mu,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=liData%real1
    mu=liData%real2
    w = liData%real3
    CalFplusprimEnd = liData%real4

    primStar = li_efold_primitive(x,alpha,mu)
    epsOneStar = li_epsilon_one(x,alpha,mu)
    potStar = li_norm_potential(x,alpha,mu)

    find_li_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_li_x_star



  function li_lnrhoend(alpha,mu,Pstar) 
    implicit none
    real(kp) :: li_lnrhoend
    real(kp), intent(in) :: alpha,mu,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = li_x_endinf(alpha,mu)
    potEnd  = li_norm_potential(xEnd,alpha,mu)
    epsOneEnd = li_epsilon_one(xEnd,alpha,mu)
    

!   Trick to return x such that rho_reh=rho_end

    x = li_x_star(alpha,mu,wrad,junk,Pstar)    
    potStar = li_norm_potential(x,alpha,mu)
    epsOneStar = li_epsilon_one(x,alpha,mu)

    
    if (.not.slowroll_validity(epsOneStar)) stop 'li_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    li_lnrhoend = lnRhoEnd

  end function li_lnrhoend

  
end module lireheat

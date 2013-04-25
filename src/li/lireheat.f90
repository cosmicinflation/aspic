!loop inflation reheating functions in the slow-roll approximations

module lireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use specialinf, only : lambert
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use lisr, only : li_epsilon_one, li_epsilon_two, li_epsilon_three
  use lisr, only : li_norm_potential
  use lisr, only : li_x_endinf, li_efold_primitive
  implicit none

  private

  public li_x_star, li_lnrhoreh_max 

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function li_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: li_x_star
    real(kp), intent(in) :: alpha,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: liData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = li_x_endinf(alpha)
    epsOneEnd = li_epsilon_one(xEnd,alpha)
    potEnd = li_norm_potential(xEnd,alpha)
    primEnd = li_efold_primitive(xEnd,alpha)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    liData%real1 = alpha
    liData%real2 = w
    liData%real3 = calF + primEnd

    if (alpha .gt. 0._kp) then
    mini = xEnd
    maxi = xEnd*1000._kp
    else
    mini = -1._kp/sqrt(2._kp) &
         /lambert(-exp(1._kp/alpha)/(sqrt(2._kp)),-1) !location of the other solution of epsilon1=1
    maxi = xEnd
    endif

    x = zbrent(find_li_x_star,mini,maxi,tolzbrent,liData)
    li_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (li_efold_primitive(x,alpha) - primEnd)
    endif

  end function li_x_star

  function find_li_x_star(x,liData)   
    implicit none
    real(kp) :: find_li_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: liData

    real(kp) :: primStar,alpha,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=liData%real1
    w = liData%real2
    CalFplusprimEnd = liData%real3

    primStar = li_efold_primitive(x,alpha)
    epsOneStar = li_epsilon_one(x,alpha)
    potStar = li_norm_potential(x,alpha)

    find_li_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_li_x_star



  function li_lnrhoreh_max(alpha,Pstar) 
    implicit none
    real(kp) :: li_lnrhoreh_max
    real(kp), intent(in) :: alpha,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = li_x_endinf(alpha)
    potEnd  = li_norm_potential(xEnd,alpha)
    epsOneEnd = li_epsilon_one(xEnd,alpha)
    

!   Trick to return x such that rho_reh=rho_end

    x = li_x_star(alpha,wrad,junk,Pstar)    
    potStar = li_norm_potential(x,alpha)
    epsOneStar = li_epsilon_one(x,alpha)

    
    if (.not.slowroll_validity(epsOneStar)) stop 'li_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    li_lnrhoreh_max = lnRhoEnd

  end function li_lnrhoreh_max

  
end module lireheat

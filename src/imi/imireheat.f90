!inverse monomial model reheating functions in the slow-roll approximations

module imireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use imisr, only : imi_epsilon_one, imi_epsilon_two, imi_epsilon_three
  use imisr, only : imi_norm_potential,imi_x_epsoneunity
  use imisr, only : imi_efold_primitive
  implicit none

  private

  public imi_x_star, imi_lnrhoreh_max
  public imi_x_rrad, imi_x_rreh

contains

!returns x given potential parameters, scalar power, wreh and
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

    mini = imi_x_epsoneunity(p)
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


!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function imi_x_rrad(p,xend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: imi_x_rrad
    real(kp), intent(in) :: p,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: imiData
    
    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = imi_epsilon_one(xEnd,p)
    potEnd = imi_norm_potential(xEnd,p)
    primEnd = imi_efold_primitive(xEnd,p)
   
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    imiData%real1 = p    
    imiData%real2 = calF + primEnd

    mini = imi_x_epsoneunity(p)
    maxi = xEnd

    x = zbrent(find_imi_x_rrad,mini,maxi,tolzbrent,imiData)
    imi_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (imi_efold_primitive(x,p) - primEnd)
    endif

  end function imi_x_rrad

  function find_imi_x_rrad(x,imiData)   
    implicit none
    real(kp) :: find_imi_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: imiData

    real(kp) :: primStar,p,CalFplusprimEnd,potStar,epsOneStar

    p=imiData%real1
    CalFplusprimEnd = imiData%real2

    primStar = imi_efold_primitive(x,p)
    epsOneStar = imi_epsilon_one(x,p)
    potStar = imi_norm_potential(x,p)

    find_imi_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)
  
  end function find_imi_x_rrad



!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function imi_x_rreh(p,xend,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: imi_x_rreh
    real(kp), intent(in) :: p,xend,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: imiData
    
    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = imi_epsilon_one(xEnd,p)
    potEnd = imi_norm_potential(xEnd,p)
    primEnd = imi_efold_primitive(xEnd,p)
   
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    imiData%real1 = p    
    imiData%real2 = calF + primEnd

    mini = imi_x_epsoneunity(p)
    maxi = xEnd

    x = zbrent(find_imi_x_rreh,mini,maxi,tolzbrent,imiData)
    imi_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (imi_efold_primitive(x,p) - primEnd)
    endif

  end function imi_x_rreh

  function find_imi_x_rreh(x,imiData)   
    implicit none
    real(kp) :: find_imi_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: imiData

    real(kp) :: primStar,p,CalFplusprimEnd,potStar

    p=imiData%real1
    CalFplusprimEnd = imiData%real2

    primStar = imi_efold_primitive(x,p)    
    potStar = imi_norm_potential(x,p)

    find_imi_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)
  
  end function find_imi_x_rreh




  function imi_lnrhoreh_max(p,xend,Pstar) 
    implicit none
    real(kp) :: imi_lnrhoreh_max
    real(kp), intent(in) :: p,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    potEnd  = imi_norm_potential(xEnd,p)
    epsOneEnd = imi_epsilon_one(xEnd,p)

!   Trick to return x such that rho_reh=rho_end

    x = imi_x_star(p,xend,wrad,junk,Pstar)    
    potStar = imi_norm_potential(x,p)
    epsOneStar = imi_epsilon_one(x,p)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'imi_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    imi_lnrhoreh_max = lnRhoEnd

  end function imi_lnrhoreh_max

  
end module imireheat

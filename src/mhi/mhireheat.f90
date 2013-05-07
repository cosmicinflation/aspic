!Mutated Hilltop inflation reheating functions in the slow-roll approximations

module mhireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use mhisr, only : mhi_epsilon_one, mhi_epsilon_two, mhi_epsilon_three
  use mhisr, only : mhi_norm_potential
  use mhisr, only : mhi_x_endinf, mhi_efold_primitive, mhi_x_trajectory
  implicit none

  private

  public mhi_x_star, mhi_lnrhoreh_max
  public mhi_x_rrad, mhi_x_rreh

contains

!returns x=phi/mu such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function mhi_x_star(mu,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: mhi_x_star
    real(kp), intent(in) :: mu,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar
    real(kp), parameter :: efoldMax = 120._kp

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: mhiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = mhi_x_endinf(mu)

    epsOneEnd = mhi_epsilon_one(xEnd,mu)
    potEnd = mhi_norm_potential(xEnd,mu)

    primEnd = mhi_efold_primitive(xEnd,mu)

    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    mhiData%real1 = mu
    mhiData%real2 = w
    mhiData%real3 = calF + primEnd

    mini = xEnd
!Assuming bfold>-120, otherwise if one uses too much big maxi values, numerical explosion
    maxi=mhi_x_trajectory(-efoldMax,xEnd,mu)


    x = zbrent(find_mhi_x_star,mini,maxi,tolzbrent,mhiData)
    mhi_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (mhi_efold_primitive(x,mu) - primEnd)
    endif

  end function mhi_x_star

  function find_mhi_x_star(x,mhiData)   
    implicit none
    real(kp) :: find_mhi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: mhiData

    real(kp) :: primStar,mu,w,CalFplusprimEnd,potStar,epsOneStar

    mu = mhiData%real1
    w = mhiData%real2
    CalFplusprimEnd = mhiData%real3

    primStar = mhi_efold_primitive(x,mu)
    epsOneStar = mhi_epsilon_one(x,mu)
    potStar = mhi_norm_potential(x,mu)

    find_mhi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_mhi_x_star


!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
 function mhi_x_rrad(mu,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: mhi_x_rrad
    real(kp), intent(in) :: mu,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar
    real(kp), parameter :: efoldMax = 120._kp

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: mhiData
    
    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif    
    
    xEnd = mhi_x_endinf(mu)

    epsOneEnd = mhi_epsilon_one(xEnd,mu)
    potEnd = mhi_norm_potential(xEnd,mu)

    primEnd = mhi_efold_primitive(xEnd,mu)

    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    mhiData%real1 = mu
    mhiData%real2 = calF + primEnd

    mini = xEnd
!Assuming bfold>-120, otherwise if one uses too much big maxi values, numerical explosion
    maxi=mhi_x_trajectory(-efoldMax,xEnd,mu)


    x = zbrent(find_mhi_x_rrad,mini,maxi,tolzbrent,mhiData)
    mhi_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (mhi_efold_primitive(x,mu) - primEnd)
    endif

  end function mhi_x_rrad

  function find_mhi_x_rrad(x,mhiData)   
    implicit none
    real(kp) :: find_mhi_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: mhiData

    real(kp) :: primStar,mu,CalFplusprimEnd,potStar,epsOneStar

    mu = mhiData%real1
    CalFplusprimEnd = mhiData%real2

    primStar = mhi_efold_primitive(x,mu)
    epsOneStar = mhi_epsilon_one(x,mu)
    potStar = mhi_norm_potential(x,mu)

    find_mhi_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)
  
  end function find_mhi_x_rrad


!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
 function mhi_x_rreh(mu,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: mhi_x_rreh
    real(kp), intent(in) :: mu,lnRreh
    real(kp), intent(out), optional :: bfoldstar
    real(kp), parameter :: efoldMax = 120._kp

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: mhiData
    
    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif    
    
    xEnd = mhi_x_endinf(mu)

    epsOneEnd = mhi_epsilon_one(xEnd,mu)
    potEnd = mhi_norm_potential(xEnd,mu)

    primEnd = mhi_efold_primitive(xEnd,mu)

    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    mhiData%real1 = mu
    mhiData%real2 = calF + primEnd

    mini = xEnd
!Assuming bfold>-120, otherwise if one uses too much big maxi values, numerical explosion
    maxi=mhi_x_trajectory(-efoldMax,xEnd,mu)


    x = zbrent(find_mhi_x_rreh,mini,maxi,tolzbrent,mhiData)
    mhi_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (mhi_efold_primitive(x,mu) - primEnd)
    endif

  end function mhi_x_rreh

  function find_mhi_x_rreh(x,mhiData)   
    implicit none
    real(kp) :: find_mhi_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: mhiData

    real(kp) :: primStar,mu,CalFplusprimEnd,potStar

    mu = mhiData%real1
    CalFplusprimEnd = mhiData%real2

    primStar = mhi_efold_primitive(x,mu)
    potStar = mhi_norm_potential(x,mu)

    find_mhi_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)
  
  end function find_mhi_x_rreh



  function mhi_lnrhoreh_max(mu,Pstar) 
    implicit none
    real(kp) :: mhi_lnrhoreh_max
    real(kp), intent(in) :: mu,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = mhi_x_endinf(mu)

    potEnd  = mhi_norm_potential(xEnd,mu)

    epsOneEnd = mhi_epsilon_one(xEnd,mu)

!   Trick to return x such that rho_reh=rho_end

    x = mhi_x_star(mu,wrad,junk,Pstar)  

    potStar = mhi_norm_potential(x,mu)
    epsOneStar = mhi_epsilon_one(x,mu)

    
    if (.not.slowroll_validity(epsOneStar)) stop 'mhi_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    mhi_lnrhoreh_max = lnRhoEnd

  end function mhi_lnrhoreh_max

  
end module mhireheat

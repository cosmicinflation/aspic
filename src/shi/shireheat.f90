!smeared Higgs inflation reheating functions in the slow-roll approximations

module shireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use shisr, only : shi_epsilon_one, shi_epsilon_two, shi_epsilon_three
  use shisr, only : shi_norm_potential
  use shisr, only : shi_x_endinf, shi_efold_primitive
  implicit none

  private

  public shi_x_star, shi_lnrhoreh_max
  public shi_x_rrad, shi_x_rreh

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function shi_x_star(alpha,phi0,xend,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: shi_x_star
    real(kp), intent(in) :: alpha,phi0,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: shiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif    

    epsOneEnd = shi_epsilon_one(xEnd,alpha,phi0)
    potEnd = shi_norm_potential(xEnd,alpha,phi0)

    primEnd = shi_efold_primitive(xEnd,alpha,phi0)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    shiData%real1 = alpha
    shiData%real2 = phi0
    shiData%real3 = w
    shiData%real4 = calF + primEnd

    mini = tolkp
    maxi = shi_x_endinf(alpha,phi0)*(1._kp-epsilon(1._kp))

    x = zbrent(find_shi_x_star,mini,maxi,tolzbrent,shiData)
    shi_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (shi_efold_primitive(x,alpha,phi0) - primEnd)
    endif

  end function shi_x_star

  function find_shi_x_star(x,shiData)   
    implicit none
    real(kp) :: find_shi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: shiData

    real(kp) :: primStar,alpha,phi0,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=shiData%real1
    phi0=shiData%real2
    w = shiData%real3
    CalFplusprimEnd = shiData%real4

    primStar = shi_efold_primitive(x,alpha,phi0)
    epsOneStar = shi_epsilon_one(x,alpha,phi0)
    potStar = shi_norm_potential(x,alpha,phi0)

    find_shi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_shi_x_star


!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function shi_x_rrad(alpha,phi0,xend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: shi_x_rrad
    real(kp), intent(in) :: alpha,phi0,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: shiData
    
    if (lnRRad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif    

    epsOneEnd = shi_epsilon_one(xEnd,alpha,phi0)
    potEnd = shi_norm_potential(xEnd,alpha,phi0)

    primEnd = shi_efold_primitive(xEnd,alpha,phi0)
   
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    shiData%real1 = alpha
    shiData%real2 = phi0
    shiData%real3 = calF + primEnd

    mini = tolkp
    maxi = shi_x_endinf(alpha,phi0)*(1._kp-epsilon(1._kp))

    x = zbrent(find_shi_x_rrad,mini,maxi,tolzbrent,shiData)
    shi_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (shi_efold_primitive(x,alpha,phi0) - primEnd)
    endif

  end function shi_x_rrad

  function find_shi_x_rrad(x,shiData)   
    implicit none
    real(kp) :: find_shi_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: shiData

    real(kp) :: primStar,alpha,phi0,CalFplusprimEnd,potStar,epsOneStar

    alpha=shiData%real1
    phi0=shiData%real2
    CalFplusprimEnd = shiData%real3

    primStar = shi_efold_primitive(x,alpha,phi0)
    epsOneStar = shi_epsilon_one(x,alpha,phi0)
    potStar = shi_norm_potential(x,alpha,phi0)

    find_shi_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)
  
  end function find_shi_x_rrad


!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function shi_x_rreh(alpha,phi0,xend,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: shi_x_rreh
    real(kp), intent(in) :: alpha,phi0,xend,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: shiData
    
    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif    

    epsOneEnd = shi_epsilon_one(xEnd,alpha,phi0)
    potEnd = shi_norm_potential(xEnd,alpha,phi0)

    primEnd = shi_efold_primitive(xEnd,alpha,phi0)
   
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    shiData%real1 = alpha
    shiData%real2 = phi0
    shiData%real3 = calF + primEnd

    mini = tolkp
    maxi = shi_x_endinf(alpha,phi0)*(1._kp-epsilon(1._kp))

    x = zbrent(find_shi_x_rreh,mini,maxi,tolzbrent,shiData)
    shi_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (shi_efold_primitive(x,alpha,phi0) - primEnd)
    endif

  end function shi_x_rreh

  function find_shi_x_rreh(x,shiData)   
    implicit none
    real(kp) :: find_shi_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: shiData

    real(kp) :: primStar,alpha,phi0,CalFplusprimEnd,potStar

    alpha=shiData%real1
    phi0=shiData%real2
    CalFplusprimEnd = shiData%real3

    primStar = shi_efold_primitive(x,alpha,phi0)
    potStar = shi_norm_potential(x,alpha,phi0)

    find_shi_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)
  
  end function find_shi_x_rreh



  function shi_lnrhoreh_max(alpha,phi0,xend,Pstar) 
    implicit none
    real(kp) :: shi_lnrhoreh_max
    real(kp), intent(in) :: alpha,phi0,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd    

    potEnd  = shi_norm_potential(xEnd,alpha,phi0)

    epsOneEnd = shi_epsilon_one(xEnd,alpha,phi0)


!   Trick to return x such that rho_reh=rho_end

    x = shi_x_star(alpha,phi0,xend,wrad,junk,Pstar)  

 
    potStar = shi_norm_potential(x,alpha,phi0)
    epsOneStar = shi_epsilon_one(x,alpha,phi0)

    
    if (.not.slowroll_validity(epsOneStar)) stop 'shi_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    shi_lnrhoreh_max = lnRhoEnd

  end function shi_lnrhoreh_max

  
end module shireheat

!double exponential inflation reheating functions in the slow-roll approximations

module deireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use deisr, only : dei_epsilon_one, dei_epsilon_two, dei_epsilon_three
  use deisr, only : dei_norm_potential
  use deisr, only : dei_x_endinf, dei_efold_primitive
  implicit none

  private

  public dei_x_star, dei_lnrhoreh_max
  public dei_x_rrad, dei_x_rreh

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function dei_x_star(beta,phi0,xend,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: dei_x_star
    real(kp), intent(in) :: beta,phi0,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: deiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif    

    epsOneEnd = dei_epsilon_one(xEnd,beta,phi0)
    potEnd = dei_norm_potential(xEnd,beta,phi0)

    primEnd = dei_efold_primitive(xEnd,beta,phi0)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    deiData%real1 = beta
    deiData%real2 = phi0
    deiData%real3 = w
    deiData%real4 = calF + primEnd

    mini = tolkp
    maxi = dei_x_endinf(beta,phi0)*(1._kp-epsilon(1._kp))

    x = zbrent(find_dei_x_star,mini,maxi,tolzbrent,deiData)
    dei_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (dei_efold_primitive(x,beta,phi0) - primEnd)
    endif

  end function dei_x_star

  function find_dei_x_star(x,deiData)   
    implicit none
    real(kp) :: find_dei_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: deiData

    real(kp) :: primStar,beta,phi0,w,CalFplusprimEnd,potStar,epsOneStar

    beta=deiData%real1
    phi0=deiData%real2
    w = deiData%real3
    CalFplusprimEnd = deiData%real4

    primStar = dei_efold_primitive(x,beta,phi0)
    epsOneStar = dei_epsilon_one(x,beta,phi0)
    potStar = dei_norm_potential(x,beta,phi0)

    find_dei_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_dei_x_star


!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function dei_x_rrad(beta,phi0,xend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: dei_x_rrad
    real(kp), intent(in) :: beta,phi0,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: deiData
    
    if (lnRRad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif    

    epsOneEnd = dei_epsilon_one(xEnd,beta,phi0)
    potEnd = dei_norm_potential(xEnd,beta,phi0)

    primEnd = dei_efold_primitive(xEnd,beta,phi0)
   
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    deiData%real1 = beta
    deiData%real2 = phi0
    deiData%real3 = calF + primEnd

    mini = tolkp
    maxi = dei_x_endinf(beta,phi0)*(1._kp-epsilon(1._kp))

    x = zbrent(find_dei_x_rrad,mini,maxi,tolzbrent,deiData)
    dei_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (dei_efold_primitive(x,beta,phi0) - primEnd)
    endif

  end function dei_x_rrad

  function find_dei_x_rrad(x,deiData)   
    implicit none
    real(kp) :: find_dei_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: deiData

    real(kp) :: primStar,beta,phi0,CalFplusprimEnd,potStar,epsOneStar

    beta=deiData%real1
    phi0=deiData%real2
    CalFplusprimEnd = deiData%real3

    primStar = dei_efold_primitive(x,beta,phi0)
    epsOneStar = dei_epsilon_one(x,beta,phi0)
    potStar = dei_norm_potential(x,beta,phi0)

    find_dei_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)
  
  end function find_dei_x_rrad


!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function dei_x_rreh(beta,phi0,xend,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: dei_x_rreh
    real(kp), intent(in) :: beta,phi0,xend,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: deiData
    
    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif    

    epsOneEnd = dei_epsilon_one(xEnd,beta,phi0)
    potEnd = dei_norm_potential(xEnd,beta,phi0)

    primEnd = dei_efold_primitive(xEnd,beta,phi0)
   
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    deiData%real1 = beta
    deiData%real2 = phi0
    deiData%real3 = calF + primEnd

    mini = tolkp
    maxi = dei_x_endinf(beta,phi0)*(1._kp-epsilon(1._kp))

    x = zbrent(find_dei_x_rreh,mini,maxi,tolzbrent,deiData)
    dei_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (dei_efold_primitive(x,beta,phi0) - primEnd)
    endif

  end function dei_x_rreh

  function find_dei_x_rreh(x,deiData)   
    implicit none
    real(kp) :: find_dei_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: deiData

    real(kp) :: primStar,beta,phi0,CalFplusprimEnd,potStar

    beta=deiData%real1
    phi0=deiData%real2
    CalFplusprimEnd = deiData%real3

    primStar = dei_efold_primitive(x,beta,phi0)
    potStar = dei_norm_potential(x,beta,phi0)

    find_dei_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)
  
  end function find_dei_x_rreh



  function dei_lnrhoreh_max(beta,phi0,xend,Pstar) 
    implicit none
    real(kp) :: dei_lnrhoreh_max
    real(kp), intent(in) :: beta,phi0,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd    

    potEnd  = dei_norm_potential(xEnd,beta,phi0)

    epsOneEnd = dei_epsilon_one(xEnd,beta,phi0)


!   Trick to return x such that rho_reh=rho_end

    x = dei_x_star(beta,phi0,xend,wrad,junk,Pstar)  

 
    potStar = dei_norm_potential(x,beta,phi0)
    epsOneStar = dei_epsilon_one(x,beta,phi0)

    
    if (.not.slowroll_validity(epsOneStar)) stop 'dei_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    dei_lnrhoreh_max = lnRhoEnd

  end function dei_lnrhoreh_max

  
end module deireheat

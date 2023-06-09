!radiatively corrected quartic inflation reheating functions in the slow-roll approximations

module hf1ireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
   use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use hf1isr, only : hf1i_epsilon_one, hf1i_epsilon_two, hf1i_epsilon_three
  use hf1isr, only : hf1i_norm_potential
  use hf1isr, only : hf1i_x_endinf, hf1i_efold_primitive
  implicit none

  private

  public hf1i_x_star, hf1i_lnrhoreh_max
  public hf1i_x_rrad, hf1i_x_rreh

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function hf1i_x_star(A1,xend,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: hf1i_x_star
    real(kp), intent(in) :: A1,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: hf1iData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = hf1i_epsilon_one(xEnd,A1)
    potEnd = hf1i_norm_potential(xEnd,A1)
    primEnd = hf1i_efold_primitive(xEnd,A1)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    hf1iData%real1 = A1    
    hf1iData%real2 = w
    hf1iData%real3 = calF + primEnd

    mini = xEnd
    maxi = 1._kp/epsilon(1._kp)

    x = zbrent(find_hf1i_x_star,mini,maxi,tolzbrent,hf1iData)
    hf1i_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (hf1i_efold_primitive(x,A1) - primEnd)
    endif

  end function hf1i_x_star

  function find_hf1i_x_star(x,hf1iData)   
    implicit none
    real(kp) :: find_hf1i_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: hf1iData

    real(kp) :: primStar,A1,w,CalFplusprimEnd,potStar,epsOneStar

    A1=hf1iData%real1
    w = hf1iData%real2
    CalFplusprimEnd = hf1iData%real3

    primStar = hf1i_efold_primitive(x,A1)
    epsOneStar = hf1i_epsilon_one(x,A1)
    potStar = hf1i_norm_potential(x,A1)

    find_hf1i_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_hf1i_x_star

!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function hf1i_x_rrad(A1,xend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: hf1i_x_rrad
    real(kp), intent(in) :: A1,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: hf1iData
    

    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = hf1i_epsilon_one(xEnd,A1)
    potEnd = hf1i_norm_potential(xEnd,A1)
    primEnd = hf1i_efold_primitive(xEnd,A1)
   
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    hf1iData%real1 = A1    
    hf1iData%real2 = calF + primEnd

    mini = xEnd
    maxi = 1._kp/epsilon(1._kp)

    x = zbrent(find_hf1i_x_rrad,mini,maxi,tolzbrent,hf1iData)
    hf1i_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (hf1i_efold_primitive(x,A1) - primEnd)
    endif

  end function hf1i_x_rrad

  function find_hf1i_x_rrad(x,hf1iData)   
    implicit none
    real(kp) :: find_hf1i_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: hf1iData

    real(kp) :: primStar,A1,CalFplusprimEnd,potStar,epsOneStar

    A1=hf1iData%real1
    CalFplusprimEnd = hf1iData%real2

    primStar = hf1i_efold_primitive(x,A1)
    epsOneStar = hf1i_epsilon_one(x,A1)
    potStar = hf1i_norm_potential(x,A1)

    find_hf1i_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)
  
  end function find_hf1i_x_rrad


  
!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function hf1i_x_rreh(A1,xend,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: hf1i_x_rreh
    real(kp), intent(in) :: A1,xend,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: hf1iData
    

    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = hf1i_epsilon_one(xEnd,A1)
    potEnd = hf1i_norm_potential(xEnd,A1)
    primEnd = hf1i_efold_primitive(xEnd,A1)
   
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    hf1iData%real1 = A1    
    hf1iData%real2 = calF + primEnd

    mini = xEnd
    maxi = 1._kp/epsilon(1._kp)

    x = zbrent(find_hf1i_x_rreh,mini,maxi,tolzbrent,hf1iData)
    hf1i_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (hf1i_efold_primitive(x,A1) - primEnd)
    endif

  end function hf1i_x_rreh

  function find_hf1i_x_rreh(x,hf1iData)   
    implicit none
    real(kp) :: find_hf1i_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: hf1iData

    real(kp) :: primStar,A1,CalFplusprimEnd,potStar

    A1=hf1iData%real1
    CalFplusprimEnd = hf1iData%real2

    primStar = hf1i_efold_primitive(x,A1)    
    potStar = hf1i_norm_potential(x,A1)

    find_hf1i_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)
  
  end function find_hf1i_x_rreh





  function hf1i_lnrhoreh_max(A1,xend,Pstar) 
    implicit none
    real(kp) :: hf1i_lnrhoreh_max
    real(kp), intent(in) :: A1,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd
    
    potEnd  = hf1i_norm_potential(xEnd,A1)
    epsOneEnd = hf1i_epsilon_one(xEnd,A1)

!   Trick to return x such that rho_reh=rho_end

    x = hf1i_x_star(A1,xend,wrad,junk,Pstar)    
    potStar = hf1i_norm_potential(x,A1)
    epsOneStar = hf1i_epsilon_one(x,A1)

    
    if (.not.slowroll_validity(epsOneStar)) stop 'hf1i_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    hf1i_lnrhoreh_max = lnRhoEnd

  end function hf1i_lnrhoreh_max

  
end module hf1ireheat

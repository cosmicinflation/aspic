!constant ns A inflation reheating functions in the slow-roll approximations

module cnaireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use cnaisr, only : cnai_epsilon_one, cnai_epsilon_two, cnai_epsilon_three
  use cnaisr, only : cnai_norm_potential
  use cnaisr, only : cnai_x_endinf, cnai_efold_primitive
  implicit none

  private

  public cnai_x_star, cnai_lnrhoreh_max
  public cnai_x_rrad, cnai_x_rreh

contains

!returns x such potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function cnai_x_star(alpha,xend,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: cnai_x_star
    real(kp), intent(in) :: alpha,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: cnaiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = cnai_epsilon_one(xEnd,alpha)
    potEnd = cnai_norm_potential(xEnd,alpha)
    primEnd = cnai_efold_primitive(xEnd,alpha)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    cnaiData%real1 = alpha    
    cnaiData%real2 = w
    cnaiData%real3 = calF + primEnd

    mini = epsilon(1._kp)
    maxi = xEnd*(1._kp-epsilon(1._kp))

    x = zbrent(find_cnai_x_star,mini,maxi,tolzbrent,cnaiData)
    cnai_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (cnai_efold_primitive(x,alpha) - primEnd)
    endif

  end function cnai_x_star

  function find_cnai_x_star(x,cnaiData)   
    implicit none
    real(kp) :: find_cnai_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: cnaiData

    real(kp) :: primStar,alpha,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=cnaiData%real1
    w = cnaiData%real2
    CalFplusprimEnd = cnaiData%real3

    primStar = cnai_efold_primitive(x,alpha)
    epsOneStar = cnai_epsilon_one(x,alpha)
    potStar = cnai_norm_potential(x,alpha)

    find_cnai_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_cnai_x_star


!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function cnai_x_rrad(alpha,xend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: cnai_x_rrad
    real(kp), intent(in) :: alpha,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: cnaiData
    

    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = cnai_epsilon_one(xEnd,alpha)
    potEnd = cnai_norm_potential(xEnd,alpha)
    primEnd = cnai_efold_primitive(xEnd,alpha)
   
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    cnaiData%real1 = alpha    
    cnaiData%real2 = calF + primEnd

    mini = epsilon(1._kp)
    maxi = xEnd*(1._kp-epsilon(1._kp))

    x = zbrent(find_cnai_x_rrad,mini,maxi,tolzbrent,cnaiData)
    cnai_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (cnai_efold_primitive(x,alpha) - primEnd)
    endif

  end function cnai_x_rrad

  function find_cnai_x_rrad(x,cnaiData)   
    implicit none
    real(kp) :: find_cnai_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: cnaiData

    real(kp) :: primStar,alpha,CalFplusprimEnd,potStar,epsOneStar

    alpha=cnaiData%real1
    CalFplusprimEnd = cnaiData%real2

    primStar = cnai_efold_primitive(x,alpha)
    epsOneStar = cnai_epsilon_one(x,alpha)
    potStar = cnai_norm_potential(x,alpha)

    find_cnai_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd &
         ,epsOneStar,potStar)
  
  end function find_cnai_x_rrad


!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function cnai_x_rreh(alpha,xend,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: cnai_x_rreh
    real(kp), intent(in) :: alpha,xend,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: cnaiData
    

    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = cnai_epsilon_one(xEnd,alpha)
    potEnd = cnai_norm_potential(xEnd,alpha)
    primEnd = cnai_efold_primitive(xEnd,alpha)
   
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    cnaiData%real1 = alpha    
    cnaiData%real2 = calF + primEnd

    mini = epsilon(1._kp)
    maxi = xEnd*(1._kp-epsilon(1._kp))

    x = zbrent(find_cnai_x_rreh,mini,maxi,tolzbrent,cnaiData)
    cnai_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (cnai_efold_primitive(x,alpha) - primEnd)
    endif

  end function cnai_x_rreh


  function find_cnai_x_rreh(x,cnaiData)   
    implicit none
    real(kp) :: find_cnai_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: cnaiData

    real(kp) :: primStar,alpha,CalFplusprimEnd,potStar

    alpha=cnaiData%real1
    CalFplusprimEnd = cnaiData%real2

    primStar = cnai_efold_primitive(x,alpha)    
    potStar = cnai_norm_potential(x,alpha)

    find_cnai_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd &
         ,potStar)
  
  end function find_cnai_x_rreh



  function cnai_lnrhoreh_max(alpha,xend,Pstar) 
    implicit none
    real(kp) :: cnai_lnrhoreh_max
    real(kp), intent(in) :: alpha,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd
    
    potEnd  = cnai_norm_potential(xEnd,alpha)
    epsOneEnd = cnai_epsilon_one(xEnd,alpha)

!   Trick to return x such that rho_reh=rho_end

    x = cnai_x_star(alpha,xend,wrad,junk,Pstar)    
    potStar = cnai_norm_potential(x,alpha)
    epsOneStar = cnai_epsilon_one(x,alpha)

    
    if (.not.slowroll_validity(epsOneStar)) stop 'cnai_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    cnai_lnrhoreh_max = lnRhoEnd

  end function cnai_lnrhoreh_max

  
end module cnaireheat

!constant ns B reheating functions in the slow-roll approximations

module cnbireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use cnbisr, only : cnbi_epsilon_one, cnbi_epsilon_two, cnbi_epsilon_three
  use cnbisr, only : cnbi_norm_potential, cnbi_x_epsoneunity
  use cnbisr, only : cnbi_x_endinf, cnbi_efold_primitive, cnbi_x_trajectory
  implicit none

  private

  public cnbi_x_star, cnbi_lnrhoreh_max 
  public cnbi_x_rrad, cnbi_x_rreh

contains

!returns x=phi/alpha given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function cnbi_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: cnbi_x_star
    real(kp), intent(in) :: alpha,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd
    
    real(kp), dimension(2) :: xEps1

    type(transfert) :: cnbiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = cnbi_x_endinf(alpha)
    xEps1 = cnbi_x_epsoneunity(alpha)

    epsOneEnd = cnbi_epsilon_one(xEnd,alpha)
    potEnd = cnbi_norm_potential(xEnd,alpha)

    primEnd = cnbi_efold_primitive(xEnd,alpha)

    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    cnbiData%real1 = alpha
    cnbiData%real2 = w
    cnbiData%real3 = calF + primEnd

    mini = xEnd + epsilon(1._kp)
    maxi = xEps1(2) + epsilon(1._kp)


    x = zbrent(find_cnbi_x_star,mini,maxi,tolzbrent,cnbiData)
    cnbi_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (cnbi_efold_primitive(x,alpha) - primEnd)
    endif

  end function cnbi_x_star

  function find_cnbi_x_star(x,cnbiData)   
    implicit none
    real(kp) :: find_cnbi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: cnbiData

    real(kp) :: primStar,alpha,w,CalFplusprimEnd,potStar,epsOneStar

    alpha = cnbiData%real1
    w = cnbiData%real2
    CalFplusprimEnd = cnbiData%real3

    primStar = cnbi_efold_primitive(x,alpha)
    epsOneStar = cnbi_epsilon_one(x,alpha)
    potStar = cnbi_norm_potential(x,alpha)

    find_cnbi_x_star = find_reheat(primStar,calFplusprimEnd,w &
         ,epsOneStar,potStar)
  
  end function find_cnbi_x_star

!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function cnbi_x_rrad(alpha,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: cnbi_x_rrad
    real(kp), intent(in) :: alpha,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd
    
    real(kp), dimension(2) :: xEps1

    type(transfert) :: cnbiData
    

    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = cnbi_x_endinf(alpha)
    xEps1 = cnbi_x_epsoneunity(alpha)

    epsOneEnd = cnbi_epsilon_one(xEnd,alpha)
    potEnd = cnbi_norm_potential(xEnd,alpha)

    primEnd = cnbi_efold_primitive(xEnd,alpha)

    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    cnbiData%real1 = alpha
    cnbiData%real2 = calF + primEnd

    mini = xEnd + epsilon(1._kp)
    maxi = xEps1(2) + epsilon(1._kp)


    x = zbrent(find_cnbi_x_rrad,mini,maxi,tolzbrent,cnbiData)
    cnbi_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (cnbi_efold_primitive(x,alpha) - primEnd)
    endif

  end function cnbi_x_rrad

  function find_cnbi_x_rrad(x,cnbiData)   
    implicit none
    real(kp) :: find_cnbi_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: cnbiData

    real(kp) :: primStar,alpha,w,CalFplusprimEnd,potStar,epsOneStar

    alpha = cnbiData%real1
    CalFplusprimEnd = cnbiData%real2

    primStar = cnbi_efold_primitive(x,alpha)
    epsOneStar = cnbi_epsilon_one(x,alpha)
    potStar = cnbi_norm_potential(x,alpha)

    find_cnbi_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd &
         ,epsOneStar,potStar)
  
  end function find_cnbi_x_rrad



!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function cnbi_x_rreh(alpha,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: cnbi_x_rreh
    real(kp), intent(in) :: alpha,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd
    
    real(kp), dimension(2) :: xEps1

    type(transfert) :: cnbiData
    

    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = cnbi_x_endinf(alpha)
    xEps1 = cnbi_x_epsoneunity(alpha)

    epsOneEnd = cnbi_epsilon_one(xEnd,alpha)
    potEnd = cnbi_norm_potential(xEnd,alpha)

    primEnd = cnbi_efold_primitive(xEnd,alpha)

    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    cnbiData%real1 = alpha
    cnbiData%real2 = calF + primEnd

    mini = xEnd + epsilon(1._kp)
    maxi = xEps1(2) + epsilon(1._kp)


    x = zbrent(find_cnbi_x_rreh,mini,maxi,tolzbrent,cnbiData)
    cnbi_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (cnbi_efold_primitive(x,alpha) - primEnd)
    endif

  end function cnbi_x_rreh

  function find_cnbi_x_rreh(x,cnbiData)   
    implicit none
    real(kp) :: find_cnbi_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: cnbiData

    real(kp) :: primStar,alpha,w,CalFplusprimEnd,potStar

    alpha = cnbiData%real1
    CalFplusprimEnd = cnbiData%real2

    primStar = cnbi_efold_primitive(x,alpha)   
    potStar = cnbi_norm_potential(x,alpha)

    find_cnbi_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd &
         ,potStar)
  
  end function find_cnbi_x_rreh



  function cnbi_lnrhoreh_max(alpha,Pstar) 
    implicit none
    real(kp) :: cnbi_lnrhoreh_max
    real(kp), intent(in) :: alpha,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = cnbi_x_endinf(alpha)

    potEnd  = cnbi_norm_potential(xEnd,alpha)

    epsOneEnd = cnbi_epsilon_one(xEnd,alpha)



!   Trick to return x such that rho_reh=rho_end

    x = cnbi_x_star(alpha,wrad,junk,Pstar)  

    potStar = cnbi_norm_potential(x,alpha)
    epsOneStar = cnbi_epsilon_one(x,alpha)

    if (.not.slowroll_validity(epsOneStar)) then
       write(*,*)'eps1 = ',epsOneStar
       stop 'cnbi_lnrhoreh_max: slow-roll violated!'
    endif

    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    cnbi_lnrhoreh_max = lnRhoEnd

  end function cnbi_lnrhoreh_max

  
end module cnbireheat

!kälher moduli inflation reheating functions in the slow-roll approximations

module kmiireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use kmiisr, only : kmii_epsilon_one, kmii_epsilon_two, kmii_epsilon_three
  use kmiisr, only : kmii_norm_potential
  use kmiisr, only : kmii_x_endinf, kmii_efold_primitive
  implicit none

  private

  public kmii_x_star, kmii_lnrhoreh_max
  public kmii_x_rrad, kmii_x_rreh

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function kmii_x_star(alpha,xend,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: kmii_x_star
    real(kp), intent(in) :: alpha,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: kmiiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = kmii_epsilon_one(xEnd,alpha)
    potEnd = kmii_norm_potential(xEnd,alpha)
    primEnd = kmii_efold_primitive(xEnd,alpha)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    kmiiData%real1 = alpha    
    kmiiData%real2 = w
    kmiiData%real3 = calF + primEnd

    !Assuming that inflation proceeds from the right to the left,
    !from initial high values of the field compared with the Planck mass
    mini = xEnd
    maxi = 1._kp/epsilon(1._kp)
  
    x = zbrent(find_kmii_x_star,mini,maxi,tolzbrent,kmiiData)
    kmii_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (kmii_efold_primitive(x,alpha) - primEnd)
    endif

  end function kmii_x_star

  function find_kmii_x_star(x,kmiiData)   
    implicit none
    real(kp) :: find_kmii_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: kmiiData

    real(kp) :: primStar,alpha,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=kmiiData%real1
    w = kmiiData%real2
    CalFplusprimEnd = kmiiData%real3

    primStar = kmii_efold_primitive(x,alpha)
    epsOneStar = kmii_epsilon_one(x,alpha)
    potStar = kmii_norm_potential(x,alpha)

    find_kmii_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_kmii_x_star


!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function kmii_x_rrad(alpha,xend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: kmii_x_rrad
    real(kp), intent(in) :: alpha,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: kmiiData
    
    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = kmii_epsilon_one(xEnd,alpha)
    potEnd = kmii_norm_potential(xEnd,alpha)
    primEnd = kmii_efold_primitive(xEnd,alpha)
   
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    kmiiData%real1 = alpha    
    kmiiData%real2 = calF + primEnd

    !Assuming that inflation proceeds from the right to the left,
    !from initial high values of the field compared with the Planck mass
    mini = xEnd
    maxi = 1._kp/epsilon(1._kp)
  
    x = zbrent(find_kmii_x_rrad,mini,maxi,tolzbrent,kmiiData)
    kmii_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (kmii_efold_primitive(x,alpha) - primEnd)
    endif

  end function kmii_x_rrad

  function find_kmii_x_rrad(x,kmiiData)   
    implicit none
    real(kp) :: find_kmii_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: kmiiData

    real(kp) :: primStar,alpha,CalFplusprimEnd,potStar,epsOneStar

    alpha=kmiiData%real1
    CalFplusprimEnd = kmiiData%real2

    primStar = kmii_efold_primitive(x,alpha)
    epsOneStar = kmii_epsilon_one(x,alpha)
    potStar = kmii_norm_potential(x,alpha)

    find_kmii_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)
  
  end function find_kmii_x_rrad


!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function kmii_x_rreh(alpha,xend,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: kmii_x_rreh
    real(kp), intent(in) :: alpha,xend,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: kmiiData
    
    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = kmii_epsilon_one(xEnd,alpha)
    potEnd = kmii_norm_potential(xEnd,alpha)
    primEnd = kmii_efold_primitive(xEnd,alpha)
   
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    kmiiData%real1 = alpha    
    kmiiData%real2 = calF + primEnd

    !Assuming that inflation proceeds from the right to the left,
    !from initial high values of the field compared with the Planck mass
    mini = xEnd
    maxi = 1._kp/epsilon(1._kp)
  
    x = zbrent(find_kmii_x_rreh,mini,maxi,tolzbrent,kmiiData)
    kmii_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (kmii_efold_primitive(x,alpha) - primEnd)
    endif

  end function kmii_x_rreh

  function find_kmii_x_rreh(x,kmiiData)   
    implicit none
    real(kp) :: find_kmii_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: kmiiData

    real(kp) :: primStar,alpha,CalFplusprimEnd,potStar

    alpha=kmiiData%real1
    CalFplusprimEnd = kmiiData%real2

    primStar = kmii_efold_primitive(x,alpha)    
    potStar = kmii_norm_potential(x,alpha)

    find_kmii_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)
  
  end function find_kmii_x_rreh



  function kmii_lnrhoreh_max(alpha,xend,Pstar) 
    implicit none
    real(kp) :: kmii_lnrhoreh_max
    real(kp), intent(in) :: alpha,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd
    
    potEnd  = kmii_norm_potential(xEnd,alpha)
    epsOneEnd = kmii_epsilon_one(xEnd,alpha)

   print*, 'xEnd=',xEnd,'potEnd=',potEnd,'epsOneEnd=',epsOneEnd

!   Trick to return x such that rho_reh=rho_end

    x = kmii_x_star(alpha,xend,wrad,junk,Pstar)    
    potStar = kmii_norm_potential(x,alpha)
    epsOneStar = kmii_epsilon_one(x,alpha)

    
    if (.not.slowroll_validity(epsOneStar)) stop 'kmii_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    kmii_lnrhoreh_max = lnRhoEnd

  end function kmii_lnrhoreh_max

  
end module kmiireheat

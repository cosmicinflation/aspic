!constant ns C model reheating functions in the slow-roll approximations

module cncireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf,ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use cncisr, only : cnci_epsilon_one, cnci_epsilon_two, cnci_epsilon_three
  use cncisr, only : cnci_norm_potential,cnci_efold_primitive,cnci_x_epsoneunity
  implicit none

  private

  public cnci_x_star, cnci_lnrhoreh_max,find_cnci_x_star
  public cnci_x_rrad, cnci_x_rreh

contains

!returns x potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the correspoding bfoldstar
  function cnci_x_star(alpha,xend,w,lnRhoReh,Pstar,bfold)    
    implicit none
    real(kp) :: cnci_x_star
    real(kp), intent(in) :: alpha,xEnd,w,lnRhoReh,Pstar
    real(kp), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: cnciData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = cnci_epsilon_one(xEnd,alpha)
    potEnd = cnci_norm_potential(xEnd,alpha)
    primEnd = cnci_efold_primitive(xEnd,alpha)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    cnciData%real1 = alpha
    cnciData%real2 = w
    cnciData%real3 = calF + primEnd

    maxi = xend*(1._kp-epsilon(1._kp))

!Value of x such that epsilon1=1 below which inflation cannot proceed
    mini = cnci_x_epsoneunity(alpha)    

    x = zbrent(find_cnci_x_star,mini,maxi,tolFind,cnciData)
    cnci_x_star = x

    if (present(bfold)) then
       bfold = -(cnci_efold_primitive(x,alpha) - primEnd)
    endif

  end function cnci_x_star

  function find_cnci_x_star(x,cnciData)   
    implicit none
    real(kp) :: find_cnci_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: cnciData

    real(kp) :: primStar,alpha,xend,w,CalFplusPrimEnd,potStar,epsOneStar

    alpha=cnciData%real1
    w = cnciData%real2
    CalFplusPrimEnd = cnciData%real3

    primStar = cnci_efold_primitive(x,alpha)
    epsOneStar = cnci_epsilon_one(x,alpha)
    potStar = cnci_norm_potential(x,alpha)

    find_cnci_x_star = find_reheat(PrimStar,calFplusPrimEnd,w,epsOneStar,potStar)

  end function find_cnci_x_star



!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function cnci_x_rrad(alpha,xend,lnRrad,Pstar,bfold)    
    implicit none
    real(kp) :: cnci_x_rrad
    real(kp), intent(in) :: alpha,xEnd,lnRrad,Pstar
    real(kp), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: cnciData
    

    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = cnci_epsilon_one(xEnd,alpha)
    potEnd = cnci_norm_potential(xEnd,alpha)
    primEnd = cnci_efold_primitive(xEnd,alpha)
   
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    cnciData%real1 = alpha
    cnciData%real2 = calF + primEnd

    maxi = xend*(1._kp-epsilon(1._kp))

!Value of x such that epsilon1=1 below which inflation cannot proceed
    mini = cnci_x_epsoneunity(alpha)    

    x = zbrent(find_cnci_x_rrad,mini,maxi,tolFind,cnciData)
    cnci_x_rrad = x

    if (present(bfold)) then
       bfold = -(cnci_efold_primitive(x,alpha) - primEnd)
    endif

  end function cnci_x_rrad

  function find_cnci_x_rrad(x,cnciData)   
    implicit none
    real(kp) :: find_cnci_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: cnciData

    real(kp) :: primStar,alpha,xend,CalFplusPrimEnd,potStar,epsOneStar

    alpha=cnciData%real1
    CalFplusPrimEnd = cnciData%real2

    primStar = cnci_efold_primitive(x,alpha)
    epsOneStar = cnci_epsilon_one(x,alpha)
    potStar = cnci_norm_potential(x,alpha)

    find_cnci_x_rrad = find_reheat_rrad(PrimStar,calFplusPrimEnd,epsOneStar,potStar)

  end function find_cnci_x_rrad



!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function cnci_x_rreh(alpha,xend,lnRreh,bfold)    
    implicit none
    real(kp) :: cnci_x_rreh
    real(kp), intent(in) :: alpha,xEnd,lnRreh
    real(kp), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: cnciData
    

    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = cnci_epsilon_one(xEnd,alpha)
    potEnd = cnci_norm_potential(xEnd,alpha)
    primEnd = cnci_efold_primitive(xEnd,alpha)
   
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    cnciData%real1 = alpha
    cnciData%real2 = calF + primEnd

    maxi = xend*(1._kp-epsilon(1._kp))

!Value of x such that epsilon1=1 below which inflation cannot proceed
    mini = cnci_x_epsoneunity(alpha)    

    x = zbrent(find_cnci_x_rreh,mini,maxi,tolFind,cnciData)
    cnci_x_rreh = x

    if (present(bfold)) then
       bfold = -(cnci_efold_primitive(x,alpha) - primEnd)
    endif

  end function cnci_x_rreh

  function find_cnci_x_rreh(x,cnciData)   
    implicit none
    real(kp) :: find_cnci_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: cnciData

    real(kp) :: primStar,alpha,xend,CalFplusPrimEnd,potStar

    alpha=cnciData%real1
    CalFplusPrimEnd = cnciData%real2

    primStar = cnci_efold_primitive(x,alpha)
    potStar = cnci_norm_potential(x,alpha)

    find_cnci_x_rreh = find_reheat_rreh(PrimStar,calFplusPrimEnd,potStar)

  end function find_cnci_x_rreh




  function cnci_lnrhoreh_max(alpha,xEnd,Pstar) 
    implicit none
    real(kp) :: cnci_lnrhoreh_max
    real(kp), intent(in) :: alpha,xEnd,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp), parameter :: wrad = 1._kp/3._kp
    real(kp), parameter :: junk= 0._kp
    real(kp) :: lnRhoEnd
         
    potEnd  = cnci_norm_potential(xEnd,alpha)
    epsOneEnd = cnci_epsilon_one(xEnd,alpha)
       
    x = cnci_x_star(alpha,xEnd,wrad,junk,Pstar)    
    potStar = cnci_norm_potential(x,alpha)
    epsOneStar = cnci_epsilon_one(x,alpha)
    
    if (.not.slowroll_validity(epsOneStar)) then
        print*,'xstar=',x,'  epsOneStar=',epsOneStar 
        stop 'cnci_lnrhoreh_max: slow-roll violated!'
    endif
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    cnci_lnrhoreh_max = lnRhoEnd

  end function cnci_lnrhoreh_max

  
  
end module cncireheat

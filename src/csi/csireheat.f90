!constant spectrum model reheating functions in the slow-roll approximations

module csireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf,ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use csisr, only : csi_epsilon_one, csi_epsilon_two, csi_epsilon_three
  use csisr, only : csi_norm_potential,csi_efold_primitive,csi_x_epsoneunity
  implicit none

  private

  public csi_x_star, csi_lnrhoreh_max,find_csi_x_star
  public csi_x_rrad, csi_x_rreh

contains

!returns x given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the correspoding bfoldstar
  function csi_x_star(alpha,xend,w,lnRhoReh,Pstar,bfold)    
    implicit none
    real(kp) :: csi_x_star
    real(kp), intent(in) :: alpha,xEnd,w,lnRhoReh,Pstar
    real(kp), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: csiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = csi_epsilon_one(xEnd,alpha)
    potEnd = csi_norm_potential(xEnd,alpha)
    primEnd = csi_efold_primitive(xEnd,alpha)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    csiData%real1 = alpha
    csiData%real2 = w
    csiData%real3 = calF + primEnd

    mini = xend
    maxi = csi_x_epsoneunity(alpha)    !Value of x such that epsilon1=1 above which inflation cannot proceed

    x = zbrent(find_csi_x_star,mini,maxi,tolFind,csiData)
    csi_x_star = x

    if (present(bfold)) then
       bfold = -(csi_efold_primitive(x,alpha) - primEnd)
    endif

  end function csi_x_star

  function find_csi_x_star(x,csiData)   
    implicit none
    real(kp) :: find_csi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: csiData

    real(kp) :: primStar,alpha,xend,w,CalFplusPrimEnd,potStar,epsOneStar

    alpha=csiData%real1
    w = csiData%real2
    CalFplusPrimEnd = csiData%real3

    primStar = csi_efold_primitive(x,alpha)
    epsOneStar = csi_epsilon_one(x,alpha)
    potStar = csi_norm_potential(x,alpha)


    find_csi_x_star = find_reheat(PrimStar,calFplusPrimEnd,w,epsOneStar,potStar)


  end function find_csi_x_star



!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function csi_x_rrad(alpha,xend,lnRrad,Pstar,bfold)    
    implicit none
    real(kp) :: csi_x_rrad
    real(kp), intent(in) :: alpha,xEnd,lnRrad,Pstar
    real(kp), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: csiData
    

    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = csi_epsilon_one(xEnd,alpha)
    potEnd = csi_norm_potential(xEnd,alpha)
    primEnd = csi_efold_primitive(xEnd,alpha)
   
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    csiData%real1 = alpha
    csiData%real2 = calF + primEnd

    mini = xend
    maxi = csi_x_epsoneunity(alpha)

    x = zbrent(find_csi_x_rrad,mini,maxi,tolFind,csiData)
    csi_x_rrad = x

    if (present(bfold)) then
       bfold = -(csi_efold_primitive(x,alpha) - primEnd)
    endif

  end function csi_x_rrad

  function find_csi_x_rrad(x,csiData)   
    implicit none
    real(kp) :: find_csi_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: csiData

    real(kp) :: primStar,alpha,xend,CalFplusPrimEnd,potStar,epsOneStar

    alpha=csiData%real1
    CalFplusPrimEnd = csiData%real2

    primStar = csi_efold_primitive(x,alpha)
    epsOneStar = csi_epsilon_one(x,alpha)
    potStar = csi_norm_potential(x,alpha)

    find_csi_x_rrad = find_reheat_rrad(PrimStar,calFplusPrimEnd,epsOneStar,potStar)

  end function find_csi_x_rrad


!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function csi_x_rreh(alpha,xend,lnRreh,bfold)    
    implicit none
    real(kp) :: csi_x_rreh
    real(kp), intent(in) :: alpha,xEnd,lnRreh
    real(kp), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: csiData
    

    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = csi_epsilon_one(xEnd,alpha)
    potEnd = csi_norm_potential(xEnd,alpha)
    primEnd = csi_efold_primitive(xEnd,alpha)
   
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    csiData%real1 = alpha
    csiData%real2 = calF + primEnd

    mini = xend
    maxi = csi_x_epsoneunity(alpha)

    x = zbrent(find_csi_x_rreh,mini,maxi,tolFind,csiData)
    csi_x_rreh = x

    if (present(bfold)) then
       bfold = -(csi_efold_primitive(x,alpha) - primEnd)
    endif

  end function csi_x_rreh

  function find_csi_x_rreh(x,csiData)   
    implicit none
    real(kp) :: find_csi_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: csiData

    real(kp) :: primStar,alpha,xend,CalFplusPrimEnd,potStar

    alpha=csiData%real1
    CalFplusPrimEnd = csiData%real2

    primStar = csi_efold_primitive(x,alpha)
    potStar = csi_norm_potential(x,alpha)

    find_csi_x_rreh = find_reheat_rreh(PrimStar,calFplusPrimEnd,potStar)

  end function find_csi_x_rreh


  function csi_lnrhoreh_max(alpha,xEnd,Pstar) 
    implicit none
    real(kp) :: csi_lnrhoreh_max
    real(kp), intent(in) :: alpha,xEnd,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp), parameter :: wrad = 1._kp/3._kp
    real(kp), parameter :: junk= 0._kp
    real(kp) :: lnRhoEnd
         
    potEnd  = csi_norm_potential(xEnd,alpha)
    epsOneEnd = csi_epsilon_one(xEnd,alpha)
       
    x = csi_x_star(alpha,xEnd,wrad,junk,Pstar)    
    potStar = csi_norm_potential(x,alpha)
    epsOneStar = csi_epsilon_one(x,alpha)
    
    if (.not.slowroll_validity(epsOneStar)) then
        print*,'xstar=',x,'  epsOneStar=',epsOneStar 
        stop 'csi_lnrhoreh_max: slow-roll violated!'
    endif
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    csi_lnrhoreh_max = lnRhoEnd

  end function csi_lnrhoreh_max

  
  
end module csireheat

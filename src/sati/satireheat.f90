!SuperConformal alpha-attractor B reheating functions in the slow-roll approximation

module satireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf, ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use satisr, only : sati_epsilon_one, sati_epsilon_two, sati_epsilon_three
  use satisr, only : sati_norm_potential, sati_efold_primitive, sati_x_endinf
  implicit none

  private

  public sati_x_star, sati_lnrhoreh_max
  public sati_x_rrad, sati_x_rreh

contains

!returns x given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the correspoding bfoldstar
  function sati_x_star(alpha,n,xend,w,lnRhoReh,Pstar,bfold)    
    implicit none
    real(kp) :: sati_x_star
    real(kp), intent(in) :: alpha,n,xend,w,lnRhoReh,Pstar
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: satiData
    
    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = sati_epsilon_one(xEnd,alpha,n)
    potEnd = sati_norm_potential(xEnd,alpha,n)
    primEnd = sati_efold_primitive(xEnd,alpha,n)

    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    satiData%real1 = alpha
    satiData%real2 = n
    satiData%real3 = w
    satiData%real4 = calF + primEnd

    mini = xend*(1._kp+epsilon(1._kp))
    maxi = sqrt(6._kp*alpha)/2._kp*log(8._kp*n*200._kp/(3._kp*alpha))

    x = zbrent(find_sati_x_star,mini,maxi,tolFind,satiData)
    sati_x_star = x

    if (present(bfold)) then
       bfold = -(sati_efold_primitive(x,alpha,n) - primEnd)
    endif

  end function sati_x_star

  function find_sati_x_star(x,satiData)   
    implicit none
    real(kp) :: find_sati_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: satiData

    real(kp) :: primStar,alpha,n,w,CalFplusPrimEnd,potStar,epsOneStar

    alpha=satiData%real1
    n=satiData%real2
    w = satiData%real3
    CalFplusPrimEnd = satiData%real4

    primStar = sati_efold_primitive(x,alpha,n)
    epsOneStar = sati_epsilon_one(x,alpha,n)
    potStar = sati_norm_potential(x,alpha,n)

    find_sati_x_star = find_reheat(PrimStar,calFplusPrimEnd,w,epsOneStar,potStar)

  end function find_sati_x_star

!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function sati_x_rrad(alpha,n,xend,lnRrad,Pstar,bfold)    
    implicit none
    real(kp) :: sati_x_rrad
    real(kp), intent(in) :: alpha,n,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: satiData
    
    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = sati_epsilon_one(xEnd,alpha,n)
    potEnd = sati_norm_potential(xEnd,alpha,n)
    primEnd = sati_efold_primitive(xEnd,alpha,n)

    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    satiData%real1 = alpha
    satiData%real2 = n
    satiData%real3 = calF + primEnd

    mini = xend*(1._kp+epsilon(1._kp))
    maxi = sqrt(6._kp*alpha)/2._kp*log(8._kp*n*200._kp/(3._kp*alpha))

    x = zbrent(find_sati_x_rrad,mini,maxi,tolFind,satiData)
    sati_x_rrad = x

    if (present(bfold)) then
       bfold = -(sati_efold_primitive(x,alpha,n) - primEnd)
    endif

  end function sati_x_rrad

  function find_sati_x_rrad(x,satiData)   
    implicit none
    real(kp) :: find_sati_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: satiData

    real(kp) :: primStar,alpha,n,CalFplusPrimEnd,potStar,epsOneStar

    alpha=satiData%real1
    n=satiData%real2
    CalFplusPrimEnd = satiData%real3

    primStar = sati_efold_primitive(x,alpha,n)
    epsOneStar = sati_epsilon_one(x,alpha,n)
    potStar = sati_norm_potential(x,alpha,n)

    find_sati_x_rrad = find_reheat_rrad(PrimStar,calFplusPrimEnd,epsOneStar,potStar)

  end function find_sati_x_rrad


!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function sati_x_rreh(alpha,n,xend,lnRreh,bfold)    
    implicit none
    real(kp) :: sati_x_rreh
    real(kp), intent(in) :: alpha,n,xend,lnRreh
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: satiData
    
    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = sati_epsilon_one(xEnd,alpha,n)
    potEnd = sati_norm_potential(xEnd,alpha,n)
    primEnd = sati_efold_primitive(xEnd,alpha,n)

    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    satiData%real1 = alpha
    satiData%real2 = n
    satiData%real3 = calF + primEnd

    mini = xend*(1._kp+epsilon(1._kp))
    maxi = sqrt(6._kp*alpha)/2._kp*log(8._kp*n*200._kp/(3._kp*alpha))

    x = zbrent(find_sati_x_rreh,mini,maxi,tolFind,satiData)
    sati_x_rreh = x

    if (present(bfold)) then
       bfold = -(sati_efold_primitive(x,alpha,n) - primEnd)
    endif

  end function sati_x_rreh

  function find_sati_x_rreh(x,satiData)   
    implicit none
    real(kp) :: find_sati_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: satiData

    real(kp) :: primStar,alpha,n,CalFplusPrimEnd,potStar

    alpha=satiData%real1
    n=satiData%real2
    CalFplusPrimEnd = satiData%real3

    primStar = sati_efold_primitive(x,alpha,n)
    potStar = sati_norm_potential(x,alpha,n)

    find_sati_x_rreh = find_reheat_rreh(PrimStar,calFplusPrimEnd,potStar)

  end function find_sati_x_rreh


  function sati_lnrhoreh_max(alpha,n,xend,Pstar) 
    implicit none
    real(kp) :: sati_lnrhoreh_max
    real(kp), intent(in) :: alpha,n,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp), parameter :: wrad = 1._kp/3._kp
    real(kp), parameter :: junk= 0._kp
    real(kp) :: lnRhoEnd
        
    potEnd  = sati_norm_potential(xEnd,alpha,n)
    epsOneEnd = sati_epsilon_one(xEnd,alpha,n)
       
    x = sati_x_star(alpha,n,xend,wrad,junk,Pstar)

    potStar = sati_norm_potential(x,alpha,n)
    epsOneStar = sati_epsilon_one(x,alpha,n)
    
    if (.not.slowroll_validity(epsOneStar)) then
        print*,'xstar=',x,'  epsOneStar=',epsOneStar 
        stop 'sati_lnrhoreh_max: slow-roll violated!'
    endif
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    sati_lnrhoreh_max = lnRhoEnd

  end function sati_lnrhoreh_max

  
  
end module satireheat

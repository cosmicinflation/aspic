!reheating functions for the radiatively corrected plateau inflation potentials
module rcpi1reheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf, ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use rcpicommon, only : rcpi_norm_potential, rcpi_efold_primitive
  use rcpicommon, only : rcpi_epsilon_one
  use rcpi1sr, only : rcpi1_check_params
  use rcpi1sr, only : rcpi1_x_endinf, rcpi1_xinimax

  implicit none

  private

  public rcpi1_x_star, rcpi1_x_rrad, rcpi1_x_rreh
  public rcpi1_lnrhoreh_max
  
contains

  
!returns x given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the correspoding bfoldstar
  function rcpi1_x_star(p,alpha,beta,w,lnRhoReh,Pstar,bfold)    
    implicit none
    real(kp) :: rcpi1_x_star
    real(kp), intent(in) :: p,alpha,beta,w,lnRhoReh,Pstar
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xEnd,potEnd
    type(transfert) :: rcpi1Data

    if (.not.rcpi1_check_params(p,alpha,beta)) stop 'rcpi1_x_star: out-of-range params!'
    
    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif

    xEnd=rcpi1_x_endinf(p,alpha,beta)
    epsOneEnd = rcpi_epsilon_one(xEnd,p,alpha,beta)
    potEnd = rcpi_norm_potential(xEnd,p,alpha,beta)
    primEnd = rcpi_efold_primitive(xEnd,p,alpha,beta)
    
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)
    
    rcpi1Data%real1 = p
    rcpi1Data%real2 = alpha
    rcpi1Data%real3 = beta
    rcpi1Data%real4 = w
    rcpi1Data%real5 = calF + primEnd

    mini = xend
    maxi = rcpi1_xinimax(p,alpha,beta)

    x = zbrent(find_rcpi1_x_star,mini,maxi,tolFind,rcpi1Data)
    rcpi1_x_star = x

    if (present(bfold)) then
       bfold = -(rcpi_efold_primitive(x,p,alpha,beta) - primEnd)
    endif

  end function rcpi1_x_star

  
  function find_rcpi1_x_star(x,rcpi1Data)   
    implicit none
    real(kp) :: find_rcpi1_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rcpi1Data

    real(kp) :: primStar,w,CalFplusPrimEnd,potStar,epsOneStar
    real(kp) :: p, alpha, beta
    
    p = rcpi1Data%real1
    alpha = rcpi1Data%real2
    beta = rcpi1Data%real3
    w = rcpi1Data%real4
    CalFplusPrimEnd = rcpi1Data%real5

    primStar = rcpi_efold_primitive(x,p,alpha,beta)
    epsOneStar = rcpi_epsilon_one(x,p,alpha,beta)
    potStar = rcpi_norm_potential(x,p,alpha,beta)

    find_rcpi1_x_star = find_reheat(PrimStar,calFplusPrimEnd,w,epsOneStar,potStar)

  end function find_rcpi1_x_star


  
!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function rcpi1_x_rrad(p,alpha,beta,lnRrad,Pstar,bfold)    
    implicit none
    real(kp) :: rcpi1_x_rrad
    real(kp), intent(in) :: p,alpha,beta,lnRrad,Pstar
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolfind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xEnd,potEnd

    type(transfert) :: rcpi1Data

    if (.not.rcpi1_check_params(p,alpha,beta)) stop 'rcpi1_x_star: out-of-range params!'
    
    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif

        
    xEnd=rcpi1_x_endinf(p,alpha,beta)
    epsOneEnd = rcpi_epsilon_one(xEnd,p,alpha,beta)
    potEnd = rcpi_norm_potential(xEnd,p,alpha,beta)
    primEnd = rcpi_efold_primitive(xEnd,p,alpha,beta)

    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    rcpi1Data%real1 = p
    rcpi1Data%real2 = alpha
    rcpi1Data%real3 = beta
    rcpi1Data%real4 = calF + primEnd

    mini = xend
    maxi = rcpi1_xinimax(p,alpha,beta)

    x = zbrent(find_rcpi1_x_rrad,mini,maxi,tolFind,rcpi1Data)
    rcpi1_x_rrad = x

    if (present(bfold)) then
       bfold = -(rcpi_efold_primitive(x,p,alpha,beta) - primEnd)
    endif

  end function rcpi1_x_rrad

  function find_rcpi1_x_rrad(x,rcpi1Data)   
    implicit none
    real(kp) :: find_rcpi1_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rcpi1Data

    real(kp) :: p,alpha,beta
    real(kp) :: primStar,CalFplusPrimEnd,potStar,epsOneStar

    p = rcpi1Data%real1
    alpha = rcpi1Data%real2
    beta = rcpi1Data%real3
    CalFplusPrimEnd = rcpi1Data%real4

    primStar = rcpi_efold_primitive(x,p,alpha,beta)
    epsOneStar = rcpi_epsilon_one(x,p,alpha,beta)
    potStar = rcpi_norm_potential(x,p,alpha,beta)

    find_rcpi1_x_rrad = find_reheat_rrad(PrimStar,calFplusPrimEnd,epsOneStar,potStar)

  end function find_rcpi1_x_rrad


!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function rcpi1_x_rreh(p,alpha,beta,lnRreh,bfold)    
    implicit none
    real(kp) :: rcpi1_x_rreh
    real(kp), intent(in) :: p,alpha,beta,lnRreh
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xEnd,potEnd

    type(transfert) :: rcpi1Data

    if (.not.rcpi1_check_params(p,alpha,beta)) stop 'rcpi1_x_star: out-of-range params!'
    
    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = rcpi1_x_endinf(p,alpha,beta)

    epsOneEnd = rcpi_epsilon_one(xEnd,p,alpha,beta)
    potEnd = rcpi_norm_potential(xEnd,p,alpha,beta)
    primEnd = rcpi_efold_primitive(xEnd,p,alpha,beta)

    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    rcpi1Data%real1 = p
    rcpi1Data%real2 = alpha
    rcpi1Data%real3 = beta
    rcpi1Data%real4 = calF + primEnd

    mini = xend
    maxi = rcpi1_xinimax(p,alpha,beta)

    x = zbrent(find_rcpi1_x_rreh,mini,maxi,tolFind,rcpi1Data)
    rcpi1_x_rreh = x

    if (present(bfold)) then
       bfold = -(rcpi_efold_primitive(x,p,alpha,beta) - primEnd)
    endif

  end function rcpi1_x_rreh

  function find_rcpi1_x_rreh(x,rcpi1Data)   
    implicit none
    real(kp) :: find_rcpi1_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rcpi1Data

    real(kp) :: p,alpha,beta
    real(kp) :: primStar,CalFplusPrimEnd,potStar

    p = rcpi1Data%real1
    alpha = rcpi1Data%real2
    beta = rcpi1Data%real3
    CalFplusPrimEnd = rcpi1Data%real4

    primStar = rcpi_efold_primitive(x,p,alpha,beta)
    potStar = rcpi_norm_potential(x,p,alpha,beta)

    find_rcpi1_x_rreh = find_reheat_rreh(PrimStar,calFplusPrimEnd,potStar)

  end function find_rcpi1_x_rreh


  function rcpi1_lnrhoreh_max(p,alpha,beta,Pstar) 
    implicit none
    real(kp) :: rcpi1_lnrhoreh_max
    real(kp), intent(in) :: p,alpha,beta,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp), parameter :: wrad = 1._kp/3._kp
    real(kp), parameter :: junk= 0._kp
    real(kp) :: lnRhoEnd
        
    xEnd = rcpi1_x_endinf(p,alpha,beta)
    potEnd  = rcpi_norm_potential(xEnd,p,alpha,beta)
    epsOneEnd = rcpi_epsilon_one(xEnd,p,alpha,beta)
       
    x = rcpi1_x_star(p,alpha,beta,wrad,junk,Pstar)

    potStar = rcpi_norm_potential(x,p,alpha,beta)
    epsOneStar = rcpi_epsilon_one(x,p,alpha,beta)
    
    if (.not.slowroll_validity(epsOneStar)) then
        print*,'xstar=',x,'  epsOneStar=',epsOneStar 
        stop 'rcpi1_lnrhoreh_max: slow-roll violated!'
    endif
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    rcpi1_lnrhoreh_max = lnRhoEnd

  end function rcpi1_lnrhoreh_max

 

end module rcpi1reheat

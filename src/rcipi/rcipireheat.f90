!reheating functions for the radiatively corrected inflection point
!inflation potentials
module rcipireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf, ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use rcipicommon, only : rcipi_norm_potential, rcipi_efold_primitive
  use rcipicommon, only : rcipi_epsilon_one
  use rcipisr, only : rcipi_check_params
  use rcipisr, only : rcipi_x_endinf, rcipi_xinimax

  implicit none

  private

  public rcipi_x_star, rcipi_x_rrad, rcipi_x_rreh
  public rcipi_lnrhoreh_max
  
contains

  
!returns x given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the correspoding bfoldstar
  function rcipi_x_star(p,alpha,beta,xend,w,lnRhoReh,Pstar,bfold)    
    implicit none
    real(kp) :: rcipi_x_star
    real(kp), intent(in) :: p,alpha,beta,xend,w,lnRhoReh,Pstar
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd
    type(transfert) :: rcipiData

    if (.not.rcipi_check_params(p,alpha,beta)) stop 'rcipi_x_star: out-of-range params!'
    
    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = rcipi_epsilon_one(xEnd,p,alpha,beta)
    potEnd = rcipi_norm_potential(xEnd,p,alpha,beta)
    primEnd = rcipi_efold_primitive(xEnd,p,alpha,beta)
    
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)
    
    rcipiData%real1 = p
    rcipiData%real2 = alpha
    rcipiData%real3 = beta
    rcipiData%real4 = w
    rcipiData%real5 = calF + primEnd

    mini = xend
    maxi = rcipi_xinimax(p,alpha,beta)

    x = zbrent(find_rcipi_x_star,mini,maxi,tolFind,rcipiData)
    rcipi_x_star = x

    if (present(bfold)) then
       bfold = -(rcipi_efold_primitive(x,p,alpha,beta) - primEnd)
    endif

  end function rcipi_x_star

  
  function find_rcipi_x_star(x,rcipiData)   
    implicit none
    real(kp) :: find_rcipi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rcipiData

    real(kp) :: primStar,w,CalFplusPrimEnd,potStar,epsOneStar
    real(kp) :: p, alpha, beta
    
    p = rcipiData%real1
    alpha = rcipiData%real2
    beta = rcipiData%real3
    w = rcipiData%real4
    CalFplusPrimEnd = rcipiData%real5

    primStar = rcipi_efold_primitive(x,p,alpha,beta)
    epsOneStar = rcipi_epsilon_one(x,p,alpha,beta)
    potStar = rcipi_norm_potential(x,p,alpha,beta)

    find_rcipi_x_star = find_reheat(PrimStar,calFplusPrimEnd,w,epsOneStar,potStar)

  end function find_rcipi_x_star


  
!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function rcipi_x_rrad(p,alpha,beta,xend,lnRrad,Pstar,bfold)    
    implicit none
    real(kp) :: rcipi_x_rrad
    real(kp), intent(in) :: p,alpha,beta,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolfind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: rcipiData

    if (.not.rcipi_check_params(p,alpha,beta)) stop 'rcipi_x_star: out-of-range params!'
    
    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif
        
    epsOneEnd = rcipi_epsilon_one(xEnd,p,alpha,beta)
    potEnd = rcipi_norm_potential(xEnd,p,alpha,beta)
    primEnd = rcipi_efold_primitive(xEnd,p,alpha,beta)

    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    rcipiData%real1 = p
    rcipiData%real2 = alpha
    rcipiData%real3 = beta
    rcipiData%real4 = calF + primEnd

    mini = xend
    maxi = rcipi_xinimax(p,alpha,beta)

    x = zbrent(find_rcipi_x_rrad,mini,maxi,tolFind,rcipiData)
    rcipi_x_rrad = x

    if (present(bfold)) then
       bfold = -(rcipi_efold_primitive(x,p,alpha,beta) - primEnd)
    endif

  end function rcipi_x_rrad

  function find_rcipi_x_rrad(x,rcipiData)   
    implicit none
    real(kp) :: find_rcipi_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rcipiData

    real(kp) :: p,alpha,beta
    real(kp) :: primStar,CalFplusPrimEnd,potStar,epsOneStar

    p = rcipiData%real1
    alpha = rcipiData%real2
    beta = rcipiData%real3
    CalFplusPrimEnd = rcipiData%real4

    primStar = rcipi_efold_primitive(x,p,alpha,beta)
    epsOneStar = rcipi_epsilon_one(x,p,alpha,beta)
    potStar = rcipi_norm_potential(x,p,alpha,beta)

    find_rcipi_x_rrad = find_reheat_rrad(PrimStar,calFplusPrimEnd,epsOneStar,potStar)

  end function find_rcipi_x_rrad


!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function rcipi_x_rreh(p,alpha,beta,xend,lnRreh,bfold)    
    implicit none
    real(kp) :: rcipi_x_rreh
    real(kp), intent(in) :: p,alpha,beta,xend,lnRreh
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: rcipiData

    if (.not.rcipi_check_params(p,alpha,beta)) stop 'rcipi_x_star: out-of-range params!'
    
    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif    

    epsOneEnd = rcipi_epsilon_one(xEnd,p,alpha,beta)
    potEnd = rcipi_norm_potential(xEnd,p,alpha,beta)
    primEnd = rcipi_efold_primitive(xEnd,p,alpha,beta)

    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    rcipiData%real1 = p
    rcipiData%real2 = alpha
    rcipiData%real3 = beta
    rcipiData%real4 = calF + primEnd

    mini = xend
    maxi = rcipi_xinimax(p,alpha,beta)

    x = zbrent(find_rcipi_x_rreh,mini,maxi,tolFind,rcipiData)
    rcipi_x_rreh = x

    if (present(bfold)) then
       bfold = -(rcipi_efold_primitive(x,p,alpha,beta) - primEnd)
    endif

  end function rcipi_x_rreh

  function find_rcipi_x_rreh(x,rcipiData)   
    implicit none
    real(kp) :: find_rcipi_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rcipiData

    real(kp) :: p,alpha,beta
    real(kp) :: primStar,CalFplusPrimEnd,potStar

    p = rcipiData%real1
    alpha = rcipiData%real2
    beta = rcipiData%real3
    CalFplusPrimEnd = rcipiData%real4

    primStar = rcipi_efold_primitive(x,p,alpha,beta)
    potStar = rcipi_norm_potential(x,p,alpha,beta)

    find_rcipi_x_rreh = find_reheat_rreh(PrimStar,calFplusPrimEnd,potStar)

  end function find_rcipi_x_rreh


  function rcipi_lnrhoreh_max(p,alpha,beta,xend,Pstar) 
    implicit none
    real(kp) :: rcipi_lnrhoreh_max
    real(kp), intent(in) :: p,alpha,beta,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp), parameter :: wrad = 1._kp/3._kp
    real(kp), parameter :: junk= 0._kp
    real(kp) :: lnRhoEnd
        
    potEnd  = rcipi_norm_potential(xEnd,p,alpha,beta)
    epsOneEnd = rcipi_epsilon_one(xEnd,p,alpha,beta)

    x = rcipi_x_star(p,alpha,beta,xend,wrad,junk,Pstar)

    potStar = rcipi_norm_potential(x,p,alpha,beta)
    epsOneStar = rcipi_epsilon_one(x,p,alpha,beta)
    
    if (.not.slowroll_validity(epsOneStar)) then
        print*,'xstar=',x,'  epsOneStar=',epsOneStar 
        stop 'rcipi_lnrhoreh_max: slow-roll violated!'
    endif
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    rcipi_lnrhoreh_max = lnRhoEnd

  end function rcipi_lnrhoreh_max

 

end module rcipireheat

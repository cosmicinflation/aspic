!Generalized Mixed large field inflation reheating functions in the
!slow-roll approximations

module gmlfireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use gmlfisr, only : gmlfi_epsilon_one, gmlfi_epsilon_two, gmlfi_epsilon_three
  use gmlfisr, only : gmlfi_norm_potential, gmlfi_x_endinf
  use gmlfisr, only :  gmlfi_efold_primitive
  implicit none

  private

  public gmlfi_x_star, gmlfi_lnrhoreh_max
  public gmlfi_x_rrad, gmlfi_x_rreh

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function gmlfi_x_star(p,q,alpha,xend,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: gmlfi_x_star
    real(kp), intent(in) :: p,q,alpha,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: gmlfiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = gmlfi_epsilon_one(xEnd,p,q,alpha)
    potEnd = gmlfi_norm_potential(xEnd,p,q,alpha)
    primEnd = gmlfi_efold_primitive(xEnd,p,q,alpha) 
    

    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    gmlfiData%real1 = alpha
    gmlfiData%real2 = p
    gmlfiData%real3 = q
    gmlfiData%real4 = xEnd
    gmlfiData%real5 = w
    gmlfiData%real6 = calF + primEnd

    mini = xend
    maxi = xend*10._kp**(6._kp)

    x = zbrent(find_gmlfi_x_star,mini,maxi,tolzbrent,gmlfiData)
    gmlfi_x_star = x  



    if (present(bfoldstar)) then
       bfoldstar = - (gmlfi_efold_primitive(x,p,q,alpha) - primEnd)
    endif

  end function gmlfi_x_star

  function find_gmlfi_x_star(x,gmlfiData)   
    implicit none
    real(kp) :: find_gmlfi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: gmlfiData

    real(kp) :: primStar,alpha,p,q,xEnd,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=gmlfiData%real1
    p=gmlfiData%real2
    q=gmlfiData%real3
    xEnd=gmlfiData%real4
    w = gmlfiData%real5
    CalFplusprimEnd = gmlfiData%real6

    primStar = gmlfi_efold_primitive(x,p,q,alpha)
    epsOneStar = gmlfi_epsilon_one(x,p,q,alpha)
    potStar = gmlfi_norm_potential(x,p,q,alpha)

    find_gmlfi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_gmlfi_x_star


!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function gmlfi_x_rrad(p,q,alpha,xend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: gmlfi_x_rrad
    real(kp), intent(in) :: p,q,alpha,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: gmlfiData
    

    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif
        
    epsOneEnd = gmlfi_epsilon_one(xEnd,p,q,alpha)
    potEnd = gmlfi_norm_potential(xEnd,p,q,alpha)
    primEnd = gmlfi_efold_primitive(xEnd,p,q,alpha) 
    
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    gmlfiData%real1 = alpha
    gmlfiData%real2 = p
    gmlfiData%real3 = q
    gmlfiData%real4 = xEnd
    gmlfiData%real5 = calF + primEnd

    mini = xend
    maxi = xend*10._kp**(6._kp)

    x = zbrent(find_gmlfi_x_rrad,mini,maxi,tolzbrent,gmlfiData)
    gmlfi_x_rrad = x  

    if (present(bfoldstar)) then
       bfoldstar = - (gmlfi_efold_primitive(x,p,q,alpha) - primEnd)
    endif

  end function gmlfi_x_rrad

  function find_gmlfi_x_rrad(x,gmlfiData)   
    implicit none
    real(kp) :: find_gmlfi_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: gmlfiData

    real(kp) :: primStar,alpha,p,q,xEnd,CalFplusprimEnd,potStar,epsOneStar

    alpha=gmlfiData%real1
    p=gmlfiData%real2
    q=gmlfiData%real3
    xEnd=gmlfiData%real4
    CalFplusprimEnd = gmlfiData%real5

    primStar = gmlfi_efold_primitive(x,p,q,alpha)
    epsOneStar = gmlfi_epsilon_one(x,p,q,alpha)
    potStar = gmlfi_norm_potential(x,p,q,alpha)

    find_gmlfi_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)
  
  end function find_gmlfi_x_rrad



  
!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function gmlfi_x_rreh(p,q,alpha,xend,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: gmlfi_x_rreh
    real(kp), intent(in) :: p,q,alpha,xend,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: gmlfiData
    

    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif
        
    epsOneEnd = gmlfi_epsilon_one(xEnd,p,q,alpha)
    potEnd = gmlfi_norm_potential(xEnd,p,q,alpha)
    primEnd = gmlfi_efold_primitive(xEnd,p,q,alpha) 
    
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    gmlfiData%real1 = alpha
    gmlfiData%real2 = p
    gmlfiData%real3 = q
    gmlfiData%real4 = xEnd
    gmlfiData%real5 = calF + primEnd

    mini = xend
    maxi = xend*10._kp**(6._kp)

    x = zbrent(find_gmlfi_x_rreh,mini,maxi,tolzbrent,gmlfiData)
    gmlfi_x_rreh = x  

    if (present(bfoldstar)) then
       bfoldstar = - (gmlfi_efold_primitive(x,p,q,alpha) - primEnd)
    endif

  end function gmlfi_x_rreh

  function find_gmlfi_x_rreh(x,gmlfiData)   
    implicit none
    real(kp) :: find_gmlfi_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: gmlfiData

    real(kp) :: primStar,alpha,p,q,xEnd,CalFplusprimEnd,potStar

    alpha=gmlfiData%real1
    p=gmlfiData%real2
    q=gmlfiData%real3
    xEnd=gmlfiData%real4
    CalFplusprimEnd = gmlfiData%real5

    primStar = gmlfi_efold_primitive(x,p,q,alpha)
    potStar = gmlfi_norm_potential(x,p,q,alpha)

    find_gmlfi_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)
  
  end function find_gmlfi_x_rreh





  function gmlfi_lnrhoreh_max(p,q,alpha,xend,Pstar) 
    implicit none
    real(kp) :: gmlfi_lnrhoreh_max
    real(kp), intent(in) :: p,q,alpha,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd
    
    potEnd  = gmlfi_norm_potential(xEnd,p,q,alpha)
    epsOneEnd = gmlfi_epsilon_one(xEnd,p,q,alpha)

!   Trick to return x such that rho_reh=rho_end

    x = gmlfi_x_star(p,q,alpha,xend,wrad,junk,Pstar)    
    potStar = gmlfi_norm_potential(x,p,q,alpha)
    epsOneStar = gmlfi_epsilon_one(x,p,q,alpha)


    
    if (.not.slowroll_validity(epsOneStar)) stop 'gmlfi_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    gmlfi_lnrhoreh_max = lnRhoEnd

  end function gmlfi_lnrhoreh_max

  
end module gmlfireheat

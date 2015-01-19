!R+R^2+alpha R^3 inflation reheating functions in the slow-roll approximations

module ccsi1reheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use ccsi1sr, only : ccsi1_epsilon_one, ccsi1_epsilon_two, ccsi1_epsilon_three
  use ccsi1sr, only : ccsi1_norm_potential, ccsi1_x_endinf, ccsi1_efold_primitive
  use ccsi1sr, only : ccsi1_numacc_xinimax

  use ccsicommon, only : ccsiBig

  implicit none

  private
  

  public ccsi1_x_star, ccsi1_lnrhoreh_max
  public ccsi1_x_rrad, ccsi1_x_rreh

contains

!returns x =phi/Mp * sqrt(2/3)
  function ccsi1_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: ccsi1_x_star
    real(kp), intent(in) :: alpha,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=epsilon(1._kp)
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: ccsi1Data
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = ccsi1_x_endinf(alpha)

    epsOneEnd = ccsi1_epsilon_one(xEnd,alpha)
    potEnd = ccsi1_norm_potential(xEnd,alpha)

    primEnd = ccsi1_efold_primitive(xEnd,alpha)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    ccsi1Data%real1 = alpha
    ccsi1Data%real2 = w
    ccsi1Data%real3 = calF + primEnd
    

    if (alpha.eq.0._kp) then !Higgs Inflation Model (HI)

       mini = xEnd
       maxi=ccsiBig !to avoid numerical explosion

    else
       
       mini = xEnd       
       maxi = ccsi1_numacc_xinimax(alpha)

    endif

    x = zbrent(find_ccsi1_x_star,mini,maxi,tolzbrent,ccsi1Data)

!consistency check
    if ((alpha.ne.0._kp).and.(x.le.xend)) stop 'ccsi1_x_star: out of numerical accuracy!'

    ccsi1_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (ccsi1_efold_primitive(x,alpha) - primEnd)
    endif

  end function ccsi1_x_star


  function find_ccsi1_x_star(x,ccsi1Data)   
    implicit none
    real(kp) :: find_ccsi1_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: ccsi1Data

    real(kp) :: primStar,alpha,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=ccsi1Data%real1
    w = ccsi1Data%real2
    CalFplusprimEnd = ccsi1Data%real3

    primStar = ccsi1_efold_primitive(x,alpha)
    epsOneStar = ccsi1_epsilon_one(x,alpha)
    potStar = ccsi1_norm_potential(x,alpha)

    find_ccsi1_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_ccsi1_x_star



!returns y given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function ccsi1_x_rrad(alpha,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: ccsi1_x_rrad
    real(kp), intent(in) :: alpha,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=epsilon(1._kp)
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: ccsi1Data
    
    
    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif

    xEnd = ccsi1_x_endinf(alpha)

    epsOneEnd = ccsi1_epsilon_one(xEnd,alpha)
    potEnd = ccsi1_norm_potential(xEnd,alpha)

    primEnd = ccsi1_efold_primitive(xEnd,alpha)
   
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    ccsi1Data%real1 = alpha
    ccsi1Data%real2 = calF + primEnd
    
    if (alpha.eq.0._kp) then !Higgs Inflation Model (HI)

       mini = xEnd
       maxi=ccsiBig !to avoid numerical explosion

    else
       
       mini = xEnd
       maxi = ccsi1_numacc_xinimax(alpha)

    endif

    x = zbrent(find_ccsi1_x_rrad,mini,maxi,tolzbrent,ccsi1Data)

!consistency check
    if ((alpha.ne.0._kp).and.(x.le.xend)) stop 'ccsi1_x_rrad: out of numerical accuracy!'

    ccsi1_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (ccsi1_efold_primitive(x,alpha) - primEnd)
    endif

  end function ccsi1_x_rrad


  function find_ccsi1_x_rrad(x,ccsi1Data)   
    implicit none
    real(kp) :: find_ccsi1_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: ccsi1Data

    real(kp) :: primStar,alpha,CalFplusprimEnd,potStar,epsOneStar

    alpha=ccsi1Data%real1
    CalFplusprimEnd = ccsi1Data%real2

    primStar = ccsi1_efold_primitive(x,alpha)
    epsOneStar = ccsi1_epsilon_one(x,alpha)
    potStar = ccsi1_norm_potential(x,alpha)

    find_ccsi1_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)
  
  end function find_ccsi1_x_rrad



!returns y given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function ccsi1_x_rreh(alpha,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: ccsi1_x_rreh
    real(kp), intent(in) :: alpha,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=epsilon(1._kp)
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: ccsi1Data
    
    
    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif

    xEnd = ccsi1_x_endinf(alpha)

    epsOneEnd = ccsi1_epsilon_one(xEnd,alpha)
    potEnd = ccsi1_norm_potential(xEnd,alpha)

    primEnd = ccsi1_efold_primitive(xEnd,alpha)
   
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    ccsi1Data%real1 = alpha
    ccsi1Data%real2 = calF + primEnd
    
    if (alpha.eq.0._kp) then !Higgs Inflation Model (HI)

       mini = xEnd
       maxi=ccsiBig !to avoid numerical explosion

    else
       
       mini = xEnd
       maxi = ccsi1_numacc_xinimax(alpha)

    endif

    x = zbrent(find_ccsi1_x_rreh,mini,maxi,tolzbrent,ccsi1Data)

!consistency check
    if ((alpha.ne.0._kp).and.(x.le.xend)) stop 'ccsi1_x_rreh: out of numerical accuracy!'

    ccsi1_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (ccsi1_efold_primitive(x,alpha) - primEnd)
    endif

  end function ccsi1_x_rreh


  function find_ccsi1_x_rreh(x,ccsi1Data)
    implicit none
    real(kp) :: find_ccsi1_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: ccsi1Data

    real(kp) :: primStar,alpha,CalFplusprimEnd,potStar

    alpha=ccsi1Data%real1
    CalFplusprimEnd = ccsi1Data%real2

    primStar = ccsi1_efold_primitive(x,alpha)
    potStar = ccsi1_norm_potential(x,alpha)

    find_ccsi1_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)
  
  end function find_ccsi1_x_rreh



  function ccsi1_lnrhoreh_max(alpha,Pstar) 
    implicit none
    real(kp) :: ccsi1_lnrhoreh_max
    real(kp), intent(in) :: alpha,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd
    
    xEnd = ccsi1_x_endinf(alpha)


    potEnd  = ccsi1_norm_potential(xEnd,alpha)

    epsOneEnd = ccsi1_epsilon_one(xEnd,alpha)

!   Trick to return y such that rho_reh=rho_end

    x= ccsi1_x_star(alpha,wrad,junk,Pstar)  
     
    potStar = ccsi1_norm_potential(x,alpha)
    epsOneStar = ccsi1_epsilon_one(x,alpha)

       
    if (.not.slowroll_validity(epsOneStar)) stop 'ccsi1_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    ccsi1_lnrhoreh_max = lnRhoEnd

  end function ccsi1_lnrhoreh_max

  
end module ccsi1reheat

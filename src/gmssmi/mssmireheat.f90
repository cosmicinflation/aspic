!MSSM inflation reheating functions in the slow-roll approximations

module mssmireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use mssmisr, only : mssmi_epsilon_one, mssmi_epsilon_two, mssmi_epsilon_three
  use mssmisr, only : mssmi_norm_potential, mssmi_x_endinf, mssmi_x_epsonemin
  use mssmisr, only :  mssmi_efold_primitive
  implicit none

  private

  public mssmi_x_star, mssmi_lnrhoreh_max
  public mssmi_x_rrad, mssmi_x_rreh

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function mssmi_x_star(phi0,xend,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: mssmi_x_star
    real(kp), intent(in) :: phi0,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: mssmiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = mssmi_epsilon_one(xEnd,phi0)
    potEnd = mssmi_norm_potential(xEnd,phi0)
    primEnd = mssmi_efold_primitive(xEnd,phi0)     

    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    mssmiData%real1 = phi0 
    mssmiData%real2 = xEnd
    mssmiData%real3 = w
    mssmiData%real4 = calF + primEnd

    mini = xend
    maxi = mssmi_x_epsonemin(phi0) -epsilon(1._kp)

    x = zbrent(find_mssmi_x_star,mini,maxi,tolzbrent,mssmiData)
    mssmi_x_star = x   


    if (present(bfoldstar)) then
       bfoldstar = - (mssmi_efold_primitive(x,phi0) - primEnd)
    endif


  end function mssmi_x_star

  function find_mssmi_x_star(x,mssmiData)   
    implicit none
    real(kp) :: find_mssmi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: mssmiData

    real(kp) :: primStar,phi0,xEnd,w,CalFplusprimEnd,potStar,epsOneStar

    phi0=mssmiData%real1
    xEnd=mssmiData%real2
    w = mssmiData%real3
    CalFplusprimEnd = mssmiData%real4

    primStar = mssmi_efold_primitive(x,phi0)
    epsOneStar = mssmi_epsilon_one(x,phi0)
    potStar = mssmi_norm_potential(x,phi0)

    find_mssmi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_mssmi_x_star


!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function mssmi_x_rrad(phi0,xend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: mssmi_x_rrad
    real(kp), intent(in) :: phi0,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: mssmiData
    

    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = mssmi_epsilon_one(xEnd,phi0)
    potEnd = mssmi_norm_potential(xEnd,phi0)
    primEnd = mssmi_efold_primitive(xEnd,phi0)     

    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    mssmiData%real1 = phi0 
    mssmiData%real2 = xEnd
    mssmiData%real3 = calF + primEnd

    mini = xend
    maxi = mssmi_x_epsonemin(phi0)-epsilon(1._kp)

    x = zbrent(find_mssmi_x_rrad,mini,maxi,tolzbrent,mssmiData)
    mssmi_x_rrad = x   


    if (present(bfoldstar)) then
       bfoldstar = - (mssmi_efold_primitive(x,phi0) - primEnd)
    endif

   
  end function mssmi_x_rrad

  function find_mssmi_x_rrad(x,mssmiData)   
    implicit none
    real(kp) :: find_mssmi_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: mssmiData

    real(kp) :: primStar,phi0,xEnd,CalFplusprimEnd,potStar,epsOneStar

    phi0=mssmiData%real1
    xEnd=mssmiData%real2
    CalFplusprimEnd = mssmiData%real3

    primStar = mssmi_efold_primitive(x,phi0)
    epsOneStar = mssmi_epsilon_one(x,phi0)
    potStar = mssmi_norm_potential(x,phi0)

    find_mssmi_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)
  
  end function find_mssmi_x_rrad


!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function mssmi_x_rreh(phi0,xend,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: mssmi_x_rreh
    real(kp), intent(in) :: phi0,xend,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: mssmiData
    

    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = mssmi_epsilon_one(xEnd,phi0)
    potEnd = mssmi_norm_potential(xEnd,phi0)
    primEnd = mssmi_efold_primitive(xEnd,phi0)     

    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    mssmiData%real1 = phi0 
    mssmiData%real2 = xEnd
    mssmiData%real3 = calF + primEnd

    mini = xend
    maxi = mssmi_x_epsonemin(phi0)-epsilon(1._kp)

    x = zbrent(find_mssmi_x_rreh,mini,maxi,tolzbrent,mssmiData)
    mssmi_x_rreh = x   


    if (present(bfoldstar)) then
       bfoldstar = - (mssmi_efold_primitive(x,phi0) - primEnd)
    endif
   

  end function mssmi_x_rreh

  function find_mssmi_x_rreh(x,mssmiData)   
    implicit none
    real(kp) :: find_mssmi_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: mssmiData

    real(kp) :: primStar,phi0,xEnd,CalFplusprimEnd,potStar

    phi0=mssmiData%real1
    xEnd=mssmiData%real2
    CalFplusprimEnd = mssmiData%real3

    primStar = mssmi_efold_primitive(x,phi0)    
    potStar = mssmi_norm_potential(x,phi0)

    find_mssmi_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)
  
  end function find_mssmi_x_rreh




  function mssmi_lnrhoreh_max(phi0,xend,Pstar) 
    implicit none
    real(kp) :: mssmi_lnrhoreh_max
    real(kp), intent(in) :: phi0,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd
    
    potEnd  = mssmi_norm_potential(xEnd,phi0)
    epsOneEnd = mssmi_epsilon_one(xEnd,phi0)

!   Trick to return x such that rho_reh=rho_end

    x = mssmi_x_star(phi0,xend,wrad,junk,Pstar)    
    potStar = mssmi_norm_potential(x,phi0)
    epsOneStar = mssmi_epsilon_one(x,phi0)
   
!    if (.not.slowroll_validity(epsOneStar)) stop 'mssmi_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    mssmi_lnrhoreh_max = lnRhoEnd

  end function mssmi_lnrhoreh_max

  
end module mssmireheat

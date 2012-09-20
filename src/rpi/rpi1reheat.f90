!R-R^p inflation reheating functions in the slow-roll approximations

module rpi1reheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use rpicommon, only : rpi_x_potmax
  use rpi1sr, only : rpi1_epsilon_one, rpi1_epsilon_two, rpi1_epsilon_three
  use rpi1sr, only : rpi1_norm_potential, rpi1_x_endinf, rpi1_efold_primitive
  implicit none

  private

  public rpi1_x_star, rpi1_lnrhoend 

contains

!returns y =phi/Mp * sqrt(2/3) such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function rpi1_x_star(p,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: rpi1_x_star
    real(kp), intent(in) :: p,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=epsilon(1._kp)
    real(kp) :: mini,maxi,calF,y,yVmax
    real(kp) :: primEnd,epsOneEnd,yend,potEnd

    real(kp), parameter :: maxiMax = 100._kp

    type(transfert) :: rpi1Data
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    yEnd = rpi1_x_endinf(p)

    epsOneEnd = rpi1_epsilon_one(yEnd,p)
    potEnd = rpi1_norm_potential(yEnd,p)

    primEnd = rpi1_efold_primitive(yEnd,p)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    rpi1Data%real1 = p
    rpi1Data%real2 = w
    rpi1Data%real3 = calF + primEnd
    

    if (p.eq.1._kp) then !Higgs Inflation Model (HI)

       mini = yEnd
       maxi=maxiMax !to avoid numerical explosion

    else
       
       mini = yEnd
       yVmax = rpi_x_potmax(p)

!local maximum of the potential, numerical safety if inflation proceeds from the right to the left
       maxi = yVmax - epsilon(1._kp)
!*(1._kp-1._kp*epsilon(1._kp))

    endif

    y = zbrent(find_rpi1_x_star,mini,maxi,tolzbrent,rpi1Data)

!consistency check
    if ((p.ne.1._kp).and.(y.le.yend)) stop 'rpi1_x_star: out of numerical accuracy!'

    rpi1_x_star = y

    if (present(bfoldstar)) then
       bfoldstar = - (rpi1_efold_primitive(y,p) - primEnd)
    endif

  end function rpi1_x_star


  function find_rpi1_x_star(y,rpi1Data)   
    implicit none
    real(kp) :: find_rpi1_x_star
    real(kp), intent(in) :: y
    type(transfert), optional, intent(inout) :: rpi1Data

    real(kp) :: primStar,p,w,CalFplusprimEnd,potStar,epsOneStar

    p=rpi1Data%real1
    w = rpi1Data%real2
    CalFplusprimEnd = rpi1Data%real3

    primStar = rpi1_efold_primitive(y,p)
    epsOneStar = rpi1_epsilon_one(y,p)
    potStar = rpi1_norm_potential(y,p)

    find_rpi1_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_rpi1_x_star



  function rpi1_lnrhoend(p,Pstar) 
    implicit none
    real(kp) :: rpi1_lnrhoend
    real(kp), intent(in) :: p,Pstar

    real(kp) :: yEnd, potEnd, epsOneEnd
    real(kp) :: y, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    yEnd = rpi1_x_endinf(p)


    potEnd  = rpi1_norm_potential(yEnd,p)

    epsOneEnd = rpi1_epsilon_one(yEnd,p)

!   Trick to return y such that rho_reh=rho_end

    y = rpi1_x_star(p,wrad,junk,Pstar)  
     
    potStar = rpi1_norm_potential(y,p)
    epsOneStar = rpi1_epsilon_one(y,p)
       
    if (.not.slowroll_validity(epsOneStar)) stop 'rpi1_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    rpi1_lnrhoend = lnRhoEnd

  end function rpi1_lnrhoend

  
end module rpi1reheat

!R-R^p inflation reheating functions in the slow-roll approximations

module rpi2reheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use rpicommon, only : rpi_x_potmax
  use rpi2sr, only : rpi2_epsilon_one, rpi2_epsilon_two, rpi2_epsilon_three
  use rpi2sr, only : rpi2_norm_potential, rpi2_efold_primitive
  implicit none

  private

  public rpi2_x_star, rpi2_lnrhoend 

contains

!returns y =phi/Mp * sqrt(2/3) such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function rpi2_x_star(p,yend,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: rpi2_x_star
    real(kp), intent(in) :: p,yend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent = tolkp
    real(kp) :: mini,maxi,calF,y,yVmax
    real(kp) :: primEnd,epsOneEnd,potEnd

    real(kp), parameter :: maxiMax = 100

    type(transfert) :: rpi2Data
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
   
    epsOneEnd = rpi2_epsilon_one(yEnd,p)
    potEnd = rpi2_norm_potential(yEnd,p)

    primEnd = rpi2_efold_primitive(yEnd,p)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    rpi2Data%real1 = p
    rpi2Data%real2 = w
    rpi2Data%real3 = calF + primEnd
    

    if (p.eq.1._kp) then !Higgs Inflation Model (HI)

       mini = yEnd
       maxi=maxiMax !to avoid numerical explosion

    else

       yVmax = rpi_x_potmax(p)
       if (yend.lt.yVmax) stop 'rpi2_x_star: yend < yVmax!'
       
!local maximum of the potential, numerical safety if inflation proceeds from the left to the right
       mini = yVmax + epsilon(1._kp)
!*(1._kp+1._kp*epsilon(1._kp))
       maxi = yEnd
              
    endif


    y = zbrent(find_rpi2_x_star,mini,maxi,tolzbrent,rpi2Data)

    rpi2_x_star = y

    if (present(bfoldstar)) then
       bfoldstar = - (rpi2_efold_primitive(y,p) - primEnd)
    endif

  end function rpi2_x_star


  function find_rpi2_x_star(y,rpi2Data)   
    implicit none
    real(kp) :: find_rpi2_x_star
    real(kp), intent(in) :: y
    type(transfert), optional, intent(inout) :: rpi2Data

    real(kp) :: primStar,p,w,CalFplusprimEnd,potStar,epsOneStar

    p=rpi2Data%real1
    w = rpi2Data%real2
    CalFplusprimEnd = rpi2Data%real3

    primStar = rpi2_efold_primitive(y,p)
    epsOneStar = rpi2_epsilon_one(y,p)
    potStar = rpi2_norm_potential(y,p)

    find_rpi2_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_rpi2_x_star



  function rpi2_lnrhoend(p,yend,Pstar) 
    implicit none
    real(kp) :: rpi2_lnrhoend
    real(kp), intent(in) :: p,yend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: y, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
   
    potEnd  = rpi2_norm_potential(yEnd,p)

    epsOneEnd = rpi2_epsilon_one(yEnd,p)

!   Trick to return y such that rho_reh=rho_end

    y = rpi2_x_star(p,yend,wrad,junk,Pstar)  
     
    potStar = rpi2_norm_potential(y,p)
    epsOneStar = rpi2_epsilon_one(y,p)
   
    if (.not.slowroll_validity(epsOneStar)) stop 'rpi2_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    rpi2_lnrhoend = lnRhoEnd

  end function rpi2_lnrhoend

  
end module rpi2reheat
